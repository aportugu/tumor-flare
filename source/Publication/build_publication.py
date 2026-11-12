#!/usr/bin/env python3
"""Format canonical outputs and create the manuscript's self-contained HTML supplement.

Usage: python build_publication.py ANALYSIS_DIR PUBLICATION_DIR CODE_DIR
Only reconciled outputs supplied by the analysis pipeline are read. Source workbooks
and private reconciliation records are not copied or embedded.
"""
from __future__ import annotations
import base64
import csv
import hashlib
import html
import io
import json
import re
import sys
from collections import OrderedDict
from pathlib import Path

ANALYSIS = Path(sys.argv[1])
OUTPUT = Path(sys.argv[2])
CODE = Path(sys.argv[3]) if len(sys.argv) > 3 else ANALYSIS.parent / "Code"
HERE = Path(__file__).resolve().parent
TABLES = OUTPUT / "Tables"
TABLES.mkdir(parents=True, exist_ok=True)

def read(name: str, required: tuple[str, ...] = ()) -> list[dict[str, str]]:
    path = ANALYSIS / (name + ".csv")
    if not path.exists():
        raise FileNotFoundError(f"Missing canonical output: {path}")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or set(required) - set(reader.fieldnames):
            raise ValueError(f"{name}: missing required fields {set(required) - set(reader.fieldnames or [])}")
        return list(reader)

def number(value: str, digits: int = 1) -> str:
    if value in ("", "NA", "NaN", "Inf", "-Inf", None):
        return ""
    val = float(value)
    return f"{val:,.{digits}f}" if digits else f"{val:,.0f}"

def save(rows: list[dict[str, str]], name: str) -> Path:
    if not rows:
        raise ValueError(f"Empty display dataset: {name}")
    path = TABLES / (name + ".csv")
    columns = list(OrderedDict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path

# DISPLAY Table1 START
def format_baseline() -> list[dict[str, str]]:
    source = read("table1_baseline", ("group", "measure", "category", "n", "evaluable", "denominator", "summary"))
    source = [row for row in source if row.get("scenario", "primary") == "primary" and row["group"] != "All"]
    groups = list(OrderedDict.fromkeys(row["group"] for row in source))
    if len(groups) != 2:
        raise ValueError(f"Table1 must contain exactly two primary groups, found {groups}")
    result: OrderedDict[tuple[str, str], dict[str, str]] = OrderedDict()
    measure_order = list(OrderedDict.fromkeys(row["measure"] for row in source))
    unavailable_categories = {"", "blank", "unknown", "not applicable", "missing", "na", "unavailable"}
    for row in source:
        key = row["measure"], row["category"]
        if key not in result:
            result[key] = {"Characteristic": row["measure"], "Category": row["category"]}
        is_unavailable = row["category"].strip().lower() in unavailable_categories
        denominator = row["denominator"] if is_unavailable else (row["evaluable"] or row["denominator"])
        if row["n"] and row["denominator"] and not 0 <= float(row["n"]) <= float(row["denominator"]):
            raise ValueError("Baseline numerator exceeds total group population")
        if row["summary"]:
            value = row["summary"]
            if denominator:
                value += f"; available n={number(denominator, 0)}"
        elif row["n"] and denominator and float(denominator) > 0:
            if not 0 <= float(row["n"]) <= float(denominator):
                raise ValueError("Baseline numerator exceeds its rendered denominator")
            percentage = str(100 * float(row["n"]) / float(denominator))
            value = f"{number(row['n'], 0)}/{number(denominator, 0)} ({number(percentage)}%)"
        else:
            value = "Not evaluable (0 available)" if denominator and float(denominator) == 0 else "Not available"
        result[key][row["group"]] = value
    for row in result.values():
        if not set(groups).issubset(row):
            raise ValueError("Canonical baseline category is missing a comparison group; zero categories must be explicit")
    return sorted(result.values(), key=lambda row: (measure_order.index(row["Characteristic"]),
        row["Category"].strip().lower() in unavailable_categories,
        row["Category"].casefold()))
# DISPLAY Table1 END

# DISPLAY Table2 START
def format_clinical_course() -> list[dict[str, str]]:
    source = read("table2_clinical_course", ("section", "field", "label", "category", "n", "denominator"))
    group_key = next((key for key in ("group", "confidence_group") if key in source[0]), None)
    if group_key is None:
        raise ValueError("Table2 requires all cases / Definite-Probable / Possible group fields")
    groups = list(OrderedDict.fromkeys(row[group_key] for row in source))
    if len(groups) != 3:
        raise ValueError(f"Table2 must contain three confidence groups, found {groups}")
    denominators = {group: {row["denominator"] for row in source if row[group_key] == group} for group in groups}
    if any(len(values) != 1 for values in denominators.values()):
        raise ValueError("Table2 must use a single total case denominator within each confidence group")
    result: OrderedDict[tuple[str, str, str], dict[str, str]] = OrderedDict()
    for row in source:
        key = row["section"], row["label"], row["category"]
        if key not in result:
            result[key] = {"Domain": row["section"], "Measure": row["label"], "Category": row["category"] or "Unknown"}
        n, denom = float(row["n"]), float(row["denominator"])
        if not 0 <= n <= denom or denom <= 0:
            raise ValueError(f"Invalid Table2 denominator: {row}")
        result[key][row[group_key]] = f"{n:.0f}/{denom:.0f} ({100*n/denom:.1f}%)"
    for row in result.values():
        for group in groups:
            row.setdefault(group, f"0/{number(next(iter(denominators[group])),0)} (0.0%)")
    return list(result.values())
# DISPLAY Table2 END

# DISPLAY TableS1 START
def format_cases() -> list[dict[str, str]]:
    from public_case_contract import read_coded_cases
    return read_coded_cases(ANALYSIS / "case_coded.csv")
# DISPLAY TableS1 END

# DISPLAY TableS2 START
def format_coverage_missingness() -> list[dict[str, str]]:
    result = []
    for row in read("coverage"):
        result.append({"Domain": "Site ascertainment", "Scenario": "Primary eligibility", "Group or site": row["site"],
            "Measure": "Recorded flare ascertainment", "Total": row["eligible_courses"], "Evaluable": row["known_flare_status"],
            "Unavailable": row["unknown_flare"], "Detail": f"Master courses: {row['master_courses']}; original flare flags: {row['original_flags']}; returned records: {row['returned_records']}; returned suspected cases: {row['returned_candidates']}"})
    for row in read("missingness"):
        result.append({"Domain": "Variable availability", "Scenario": row["scenario"], "Group or site": row["group"],
            "Measure": row["field"], "Total": row["total"], "Evaluable": row["evaluable"], "Unavailable": row["unavailable"], "Detail": ""})
    return result
# DISPLAY TableS2 END

# DISPLAY TableS3 START
def format_comparisons() -> list[dict[str, str]]:
    result: list[dict[str, str]] = []
    def count_value(n: str, denom: str) -> str:
        return f"{number(n, 0)}/{number(denom, 0)} ({number(str(100 * float(n) / float(denom)))}%)" if float(denom) else "Not available"
    for name, domain, measure_key in (("response_distribution", "Best response distribution", "category"),
                                      ("response_thresholds", "Best response thresholds", "threshold")):
        combined: OrderedDict[tuple[str, str], dict[str, str]] = OrderedDict()
        for row in read(name):
            measure = {"pr": ">=PR", "vgpr": ">=VGPR", "cr": "CR/sCR"}.get(row[measure_key], row[measure_key])
            key = row["scenario"], measure
            if key not in combined:
                combined[key] = {"Domain": domain, "Scenario": row["scenario"], "Measure": measure,
                    "Flare candidate": "", "No included flare": "", "Test": "", "P value": "", "Qualification": "Best index-course response; authoritative unknowns remain unavailable"}
            numerator = row.get("responders", row.get("n", "0"))
            is_unknown = measure.lower() in ("ne/unknown", "unknown", "unavailable", "ne")
            denominator = row["total_courses"] if is_unknown else row["evaluable"]
            combined[key][row["group"]] = count_value(numerator, denominator)
        result.extend(combined.values())
    for row in read("clinical_comparisons"):
        result.append({"Domain": "Post-treatment events", "Scenario": row["scenario"], "Measure": row["label"],
            "Flare candidate": count_value(row["flare_yes"], row["flare_evaluable"]),
            "No included flare": count_value(row["comparison_yes"], row["comparison_evaluable"]),
            "Test": "Fisher exact", "P value": number(row["fisher_p"], 4), "Qualification": row["interpretation"]})
    combined = OrderedDict()
    for row in read("hospitalization"):
        key = row["scenario"]
        if key not in combined:
            combined[key] = {"Domain": "Hospital duration among admissions", "Scenario": key, "Measure": "Days, median (IQR)",
                "Flare candidate": "", "No included flare": "", "Test": "Wilcoxon rank-sum", "P value": number(row["wilcoxon_p"], 4),
                "Qualification": row["definition"] + "; " + row["test_status"]}
        combined[key][row["group"]] = f"{number(row['median'])} ({number(row['q1'])}–{number(row['q3'])}); available {number(row['admitted_evaluable'],0)}/{number(row['total_courses'],0)}"
    result.extend(combined.values())
    return result
# DISPLAY TableS3 END

# DISPLAY TableS4 START
def format_models() -> list[dict[str, str]]:
    source = read("models")
    columns = ("model_id", "scenario", "family", "label", "response", "adjustment",
               "n_courses", "n_patient_clusters", "candidates", "events", "candidate_events",
               "excluded_missing", "estimate", "conf_low", "conf_high", "p_value", "lab_bh_p",
               "bootstrap_successes", "bootstrap_requested", "bootstrap_seed", "model_status", "diagnostic", "sparse_note")
    return [{key: row.get(key, "") for key in columns} for row in source]
# DISPLAY TableS4 END

# DISPLAY TableS5 START
def format_survival() -> list[dict[str, str]]:
    result = []
    for name, domain in (("survival_models", "Time-dependent exposure"), ("landmark_models", "Landmark sensitivity")):
        for row in read(name):
            if row.get("scenario") == "original_master_endpoints":
                result.append({"Analysis": "Historical audit of superseded master endpoints", **row})
            else:
                result.append({"Analysis": domain, **row})
    return sorted(result, key=lambda row: row.get("scenario") == "original_master_endpoints")
# DISPLAY TableS5 END

TABLE_BUILDERS = OrderedDict((
    ("Table_1", format_baseline), ("Table_2", format_clinical_course),
    ("Table_S1", format_cases), ("Table_S2", format_coverage_missingness),
    ("Table_S3", format_comparisons), ("Table_S4", format_models), ("Table_S5", format_survival)))
formatted_tables = {name: builder() for name, builder in TABLE_BUILDERS.items()}
rendered_count_cells = 0
for table_name, rows in formatted_tables.items():
    for row in rows:
        for value in row.values():
            for match in re.finditer(r"\b([0-9][0-9,]*)/([0-9][0-9,]*)(?:\s+\(([0-9.]+)%\))?", value or ""):
                count, denominator = (int(match.group(i).replace(",", "")) for i in (1, 2))
                if count > denominator:
                    raise ValueError(f"Rendered count exceeds denominator in {table_name}")
                if match.group(3) is not None and not 0 <= float(match.group(3)) <= 100:
                    raise ValueError(f"Rendered percentage outside 0–100 in {table_name}")
                rendered_count_cells += 1
for name, rows in formatted_tables.items():
    save(rows, name)

CAPTIONS = OrderedDict([
    ("Table_1", ("Baseline characteristics of the comparative cohort",
       "Characteristics are described by documented suspected flare status in the eligible cohort. Observed categories use n/evaluable N and percent; missing, Blank, Unknown, and Not applicable categories use n/total group N. Continuous values show the canonical summary and available N, with unavailable counts reported over the total group population. Product, demographics, performance status, prior treatment, cytogenetics, pretreatment disease distribution, marrow involvement, secretion pattern, and baseline laboratory results precede treatment. No hypothesis tests are shown.")),
    ("Table_2", ("Management and subsequent course of site-reviewed suspected flare",
       "All site-reviewed cases, the Definite/Probable subset, and Possible cases are shown separately. Each cell uses the total number in its confidence group; unknown assessments remain a separate category. Interventions may overlap. Bispecific treatment action, symptom-resolution status, local lesion outcome, and first post-flare response during the index course are distinct measures. Responses after another systemic therapy are not attributed to the index bispecific course.")),
    ("Table_S1", ("Structured characteristics of individual site-reviewed cases",
       "Neutral case labels replace original identifiers. Every case value is a predefined code, category, or numeric elapsed value. Chronology uses days after first step-up dosing; timing codes identify recorded points, estimated points, onset intervals, and unavailable timing. Empty numeric entries mean unavailable. First and subsequent response attribution flags, symptom recovery, local lesion outcome, and best index-course systemic response remain distinct. Coded additional-treatment local assessment context identifies radiation, salvage, or CAR-T; NOT_REPORTED does not establish absence of that treatment. Clinical qualification narratives are excluded.")),
    ("Table_S2", ("Site ascertainment and variable availability",
       "Site-level flare ascertainment and the availability of baseline, response, and survival variables are shown. Courses without recorded flare ascertainment are excluded from the comparative incidence denominator. Availability is specific to each variable and analysis population; missingness is not interpreted as an absent clinical finding.")),
    ("Table_S3", ("Comparative response and post-treatment events",
       "Full best-response categories, response thresholds, CRS, ICANS, prophylactic and therapeutic tocilizumab, hospital duration, and infections within 30 days are described with their evaluable denominators. Hospital duration uses recorded total days among admitted courses and does not estimate flare-attributable duration. CRS, ICANS, supportive treatment, hospitalization, and infection are post-treatment events; these comparisons do not establish baseline risk or temporal precedence.")),
    ("Table_S4", ("Baseline laboratory and response association models",
       "All prespecified exploratory Firth logistic models and sensitivity populations are shown with course and patient-cluster counts, flare events, estimates, confidence intervals, and bootstrap diagnostics. Baseline variables are assessed independently for association with documented suspected flare. Adjusted response models use six prespecified baseline covariates. Laboratory tests include Benjamini-Hochberg adjustment. Inference follows the canonical successful-resample threshold and is suppressed when estimation criteria are not met.")),
    ("Table_S5", ("Time-dependent and landmark survival models",
       "PFS and OS models show exposure-only, corresponding complete-case, and adjusted analyses with population and event counts, confidence intervals, endpoint qualifications, and proportional-hazards diagnostics. Suspected flare becomes an exposure at recorded onset. Sensitivities examine confidence classification, source-supported alternative onset, and day-30/day-60 landmarks. Unknown authoritative endpoints are excluded from affected analyses rather than replaced with conflicting master values. The original-master-endpoint analysis is labeled as a historical audit of superseded values and is not a clinically valid sensitivity analysis.")),
    ("Figure_1", ("Construction of the comparative cohort and site-reviewed case series",
       "Panel A shows the master cohort through supported-product, first-eligible-BCMA/all-talquetamab, bridging, and recorded-ascertainment rules, followed by documented suspected flare and no documented flare groups. Panel B shows the independently reviewed case-series population, confidence assessments, and comparative eligibility exclusions. Bridging and repeat BCMA courses contribute to the descriptive case series when site reviewed but do not enter the primary comparative cohort.")),
    ("Figure_2", ("Onset and recorded clinical chronology",
       "Panel A shows the empirical cumulative distribution among cases with point onset timing. Panel B uses neutral case rows to show onset, recorded peak, and recorded symptom-resolution days after first step-up dosing. Open markers indicate estimated dates; an onset interval is drawn as an interval without an inferred midpoint. Markers beyond day 100 are placed at the right boundary with a greater-than mark; elapsed values are supplied in Table S1. A missing resolution date does not imply persistent flare. Confidence and bridging status are shown separately.")),
    ("Figure_3", ("Clinical manifestations impact and index anatomy",
       "Panels A-C describe six manifestation domains, clinical complications and flare-attributed admission, and index anatomical involvement. Domains overlap, and each uses the complete site-reviewed denominator with Yes, No, Unknown, and Not applicable categories. Unknown and Not applicable are separate states. Labels give the Yes count over the total case-series denominator. True extramedullary disease is distinguished from paraskeletal involvement. Spinal and CNS involvement are separate overlapping descriptors; spinal involvement does not by itself imply CNS disease.")),
    ("Figure_4", ("Exploratory associations with suspected flare response and survival",
       "Panel A shows six independent Firth logistic baseline models for documented suspected flare. Panel B shows adjusted best-response thresholds, and panel C shows adjusted time-dependent PFS and OS associations. Each row provides its estimate, 95% confidence interval, and actual analysis population. Baseline ORs, response ORs, and survival HRs have different outcomes and are displayed separately. Six-covariate adjustment includes age, sex, ECOG performance status, prior treatment lines, high-risk cytogenetics, and true extramedullary disease. Sparse events, missing covariates, post-treatment classification, and wide intervals limit inference; absence of a clear association does not establish equivalence.")),
    ("Figure_5", ("Proposed clinical recognition and reassessment framework",
       "The proposed framework prioritizes urgent complications, competing explanations, and integration of clinical, biochemical, and imaging evidence. Treatment decisions are individualized rather than determined by a mandatory continuation rule. Longitudinal reassessment records symptom recovery, local lesion outcome, and systemic response separately. This framework is derived from the cohort and literature and has not undergone formal consensus validation or prospective diagnostic evaluation.")),
    ("Figure_S1", ("Landmark progression-free and overall survival",
       "PFS and OS are shown at day-30 and day-60 landmarks with 95% confidence bands and numbers at risk. Follow-up restarts at the landmark. Flare exposure is defined only by onset on or before that landmark; courses developing flare later remain in the no-flare-by-landmark group for that analysis. Patients must be alive and, for PFS, progression-free at the landmark. These descriptive curves condition on landmark eligibility and do not replace the time-dependent primary models.")),
])

FIGURE_FILES = {
    "Figure_1": "Figure_1_Cohort_flow", "Figure_2": "Figure_2_Onset_and_chronology",
    "Figure_3": "Figure_3_Clinical_spectrum", "Figure_4": "Figure_4_Association_models",
    "Figure_5": "Figure_5_Recognition_framework", "Figure_S1": "Figure_S1_Landmark_survival"}
SOURCE_FILES = {
    "Table_1": ["table1_baseline.csv"], "Table_2": ["table2_clinical_course.csv"],
    "Table_S1": ["case_coded.csv", "public_case_vocabulary.json"], "Table_S2": ["coverage.csv", "missingness.csv"],
    "Table_S3": ["response_distribution.csv", "response_thresholds.csv", "clinical_comparisons.csv", "hospitalization.csv"],
    "Table_S4": ["models.csv", "bootstrap_diagnostics.csv"],
    "Table_S5": ["survival_models.csv", "landmark_models.csv"],
    "Figure_1": ["flow.csv", "incidence.csv", "case_coded.csv", "public_case_vocabulary.json"], "Figure_2": ["case_coded.csv", "public_case_vocabulary.json", "timing_summary.csv"],
    "Figure_3": ["presentation.csv", "anatomy_binary.csv"], "Figure_4": ["models.csv", "survival_models.csv"],
    "Figure_5": [], "Figure_S1": ["landmark_km_steps.csv", "landmark_risk_counts.csv"]}

# Bibliographic metadata is portable; full-text paths, abstracts, attachments,
# and native reference-manager databases are excluded from the public supplement.
reference_path = CODE.parent / "References/reference_metadata.json"
if not reference_path.exists():
    reference_path = CODE.parent / "References/reference_metadata_public.json"
reference_source = json.loads(reference_path.read_text(encoding="utf-8"))
reference_fields = ("reference_key", "first_author", "year", "authors", "title", "journal",
                    "volume", "issue", "pages_or_article", "doi", "pmid", "pmcid",
                    "pubmed_verified", "europepmc_verified", "crossref_verified", "evidence_status")
references = [{key: row.get(key, "") for key in reference_fields} for row in reference_source]
for reference in references:
    if not reference["doi"] or not re.match(r"^10\.[0-9]{4,9}/\S+$", reference["doi"], re.I):
        raise ValueError(f"Reference needs a verified DOI: {reference['reference_key']}")
    if not any(reference.get(field) for field in ("pubmed_verified", "europepmc_verified", "crossref_verified")):
        raise ValueError(f"Reference lacks bibliographic verification: {reference['reference_key']}")
reference_lookup = {row["reference_key"]: row for row in references}
(OUTPUT / "reference_metadata_public.json").write_text(json.dumps(references, indent=2, ensure_ascii=False), encoding="utf-8")
CAPTION_REFERENCES = {
    "Figure_2": ("Forsberg2024", "Nagarajan2024", "Jadot2026"),
    "Figure_3": ("Forsberg2024", "Khan2025", "Kawatsuki2026", "Ito2026", "Takakuwa2026"),
    "Figure_5": ("Forsberg2024", "Nagarajan2024", "Khan2025", "Daunov2026", "Jadot2026", "Biltibo2026", "Lee2019", "Kumar2016")}

def reference_links(keys: tuple[str, ...]) -> str:
    return "; ".join(f"<a href='#ref-{html.escape(key)}'>{html.escape(str(reference_lookup[key]['first_author']) + ' ' + str(reference_lookup[key]['year']))}</a>" for key in keys)

reference_items = []
for row in references:
    authors = row["authors"]
    author = "; ".join(authors[:3]) + ("; et al." if len(authors) > 3 else ".")
    reference_items.append(f"<li id='ref-{html.escape(row['reference_key'])}'>{html.escape(author)} {html.escape(row['title'])} <em>{html.escape(row['journal'])}</em>. {html.escape(str(row['year']))};{html.escape(str(row['volume']))}:{html.escape(str(row['pages_or_article']))}. <a href='https://doi.org/{html.escape(row['doi'])}'>doi:{html.escape(row['doi'])}</a></li>")
reference_section = "<section id='references'><h2>References</h2><p>DOIs and citation metadata were checked against bibliographic records. Direct reports support recognition of the phenotype and competing explanations; they do not validate the proposed framework.</p><ol>" + "".join(reference_items) + "</ol></section>"

def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()

def source_snippet(display: str) -> str:
    if display.startswith("Figure"):
        path = HERE / "build_figures.R"
    else:
        path = HERE / "build_publication.py"
    code = path.read_text(encoding="utf-8")
    tag = display.replace("_", "")
    match = re.search(rf"# DISPLAY {re.escape(tag)} START\n(.*?)# DISPLAY {re.escape(tag)} END", code, re.S)
    if not match:
        raise ValueError(f"Missing executable source marker for {display}")
    return match.group(1).strip()

def table_html(rows: list[dict[str, str]]) -> str:
    if "case_label" in rows[0]:
        groups = OrderedDict([
            ("Course eligibility", ("bridging", "comparative_included")),
            ("Systemic response and attribution", ("first_postflare_response", "first_response_on_index", "first_response_category_estimated", "subsequent_response", "subsequent_response_on_index", "best_index_response")),
            ("Symptom recovery and local assessment", ("symptom_status", "local_outcome", "later_local_assessment_context", "later_local_after_radiation", "later_local_after_salvage", "later_local_after_cart", "later_local_separate_from_systemic_response")),
            ("Elapsed chronology and timing status", ("onset_days", "peak_days", "resolution_days", "onset_lower", "onset_upper", "onset_estimated", "peak_estimated", "resolution_estimated", "onset_timing_type", "peak_timing_type", "resolution_timing_type", "peak_before_point_onset")),
            ("Manifestations", ("pain_flare", "palpable_enlargement", "imaging_worsening", "new_lesion", "neurologic_cns", "biochemical_or_organ", "focality", "known_or_new_site")),
            ("Anatomical involvement", ("skeletal_flare_type", "paraskeletal_involvement", "true_emd_classification", "spinal_involvement", "cns_involvement")),
            ("Clinical impact", ("organ_dysfunction", "flare_hospitalization", "icu_due_to_flare", "pathologic_fracture", "cord_compression")),
            ("Competing findings", ("overlap_crs", "overlap_icans", "overlap_infection", "alternative_fracture", "alternative_toxicity", "other_competing_explanation", "alternative_diagnosis", "true_pd_within_14d")),
            ("Management", ("bsab_action", "management_steroids", "management_radiation", "management_surgery", "management_toci_anakinra"))])
        grouped_fields = [key for fields in groups.values() for key in fields]
        if set(grouped_fields) != set(rows[0]) - {"case_label", "product", "confidence"} or len(grouped_fields) != len(set(grouped_fields)):
            raise ValueError("Case HTML groups do not cover the closed case contract exactly")
        case_tables = []
        for row in rows:
            label = f"{row['case_label']} · {row['product']} · {row['confidence']}"
            pairs = ""
            for group, fields in groups.items():
                pairs += f"<tr><th colspan='2' scope='colgroup'>{html.escape(group)}</th></tr>"
                for key in fields:
                    value = row[key]
                    display_value = {"TRUE": "Yes", "FALSE": "No"}.get(value, value or "UNAVAILABLE")
                    pairs += f"<tr><th scope='row'>{html.escape(key.replace('_', ' ').capitalize())}</th><td>{html.escape(display_value)}</td></tr>"
            case_tables.append(f"<details class='case'><summary>{html.escape(label)}</summary><table><tbody>{pairs}</tbody></table></details>")
        return "<p>Expand a neutral case label to read its coded record. The response-attribution and additional-treatment assessment fields appear first. All elapsed values use days after first step-up dosing. UNAVAILABLE numeric entries indicate missing timing; NOT_REPORTED flags do not establish absence of a treatment or assessment.</p>" + "".join(case_tables)
    keys = list(OrderedDict.fromkeys(key for row in rows for key in row))
    headings = "".join(f"<th scope='col'>{html.escape(key.replace('_', ' '))}</th>" for key in keys)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(row.get(key, '') or '')}</td>" for key in keys) + "</tr>")
    return "<div class='table-scroll'><table><thead><tr>" + headings + "</tr></thead><tbody>" + "".join(body) + "</tbody></table></div>"

def inline_png(path: Path) -> str:
    from PIL import Image
    if not path.exists():
        raise FileNotFoundError(f"Missing final figure: {path}")
    im = Image.open(path).convert("RGB")
    if im.width > 2000:
        im.thumbnail((2000, 2400), Image.Resampling.LANCZOS)
    buf = io.BytesIO()
    im.save(buf, format="PNG", optimize=True)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode("ascii")

registry = []
sections = []
for display, (title, caption) in CAPTIONS.items():
    source_names = SOURCE_FILES[display]
    for source in source_names:
        if not (ANALYSIS / source).exists():
            raise FileNotFoundError(f"{display} missing source: {source}")
    stem = display.replace("_", " ")
    if display.startswith("Table"):
        dataset = f"Tables/{display}.csv"
        body = table_html(formatted_tables[display])
        download = f"<a href='{dataset}' download>Download CSV</a>"
        module = "build_publication.py"
    else:
        filename = FIGURE_FILES[display]
        dataset = "Tables/Table_S1.csv" if display == "Figure_2" else f"Display_data/{display}.csv"
        if display == "Figure_S1": dataset = "Display_data/Figure_S1_KM.csv; Display_data/Figure_S1_Risk.csv"
        body = f"<figure><img alt='{html.escape(stem + '. ' + title)}' src='{inline_png(OUTPUT / 'Figures' / (filename + '.png'))}'></figure>"
        download = f"<a href='Figures/{filename}.pdf' download>Download PDF</a> <a href='Figures/{filename}.png' download>Download PNG</a>"
        module = "build_figures.R"
    source_hashes = "; ".join(f"{source}: {sha256(ANALYSIS / source)}" for source in source_names)
    check = "canonical numerical checks and publication display checks"
    registry.append({"display": display, "title": title, "dataset": dataset, "analysis_source": "; ".join(source_names),
                     "source_sha256": source_hashes, "publication_module": module,
                     "module_sha256": sha256(HERE / module), "verification": check, "caption": caption})
    citations = f"<p class='caption'>Relevant literature: {reference_links(CAPTION_REFERENCES[display])}.</p>" if display in CAPTION_REFERENCES else ""
    sections.append(f"<section id='{display}'><h2>{html.escape(stem + '. ' + title)}</h2><p class='caption'>{html.escape(caption)}</p>{citations}{download}{body}"
                    f"<details><summary>Executable code for {html.escape(stem)}</summary><pre><code>{html.escape(source_snippet(display))}</code></pre></details>"
                    f"<details><summary>Dataset and verification</summary><p>{html.escape('; '.join(source_names) or 'Conceptual framework; no empirical estimate')}</p>"
                    f"<p>Publication module: {html.escape(module)}. Verification: canonical analysis and display checks.</p><pre>{html.escape(source_hashes)}</pre></details></section>")
save(registry, "Display_registry")

GENERIC_R_FILES = ("run_analysis.R", "analysis_functions.R", "survival_functions.R",
                   "reconciliation_functions.R", "model_functions.R", "summary_functions.R", "ground_truth_tests.R")
code_sections = []
for filename in GENERIC_R_FILES:
    path = CODE / filename
    text = path.read_text(encoding="utf-8")
    from public_case_contract import public_generic_source
    text = public_generic_source(filename, text)
    if re.search(r"/Users/[^\s\"']+", text) or re.search(r"\b(?:CCF|TEC|TALQ|MDACC)[_-]?\d{2,}", text):
        raise ValueError(f"Private path or clinical identifier in generic source: {path.name}")
    code_sections.append(f"<details><summary>{html.escape(path.name)}</summary><pre><code>{html.escape(text)}</code></pre></details>")
if not code_sections:
    raise ValueError("No generic executable analysis source was supplied")
contract_source = (HERE / "public_case_contract.py").read_text(encoding="utf-8")
case_vocabulary = (ANALYSIS / "public_case_vocabulary.json").read_text(encoding="utf-8")
code_sections.append("<details><summary>Closed public case contract</summary><pre><code>" + html.escape(contract_source) + "</code></pre></details>")
code_sections.append("<details><summary>Finite case category vocabulary</summary><pre><code>" + html.escape(case_vocabulary) + "</code></pre></details>")
shared = """<section id='pipeline'><h2>Shared analysis pipeline</h2>
<p>The analysis starts with the master course cohort and site-reviewed electronic data collection forms. Site EDC findings take precedence for matching clinical concepts when sources conflict. Unknown authoritative values remain unknown; unresolved endpoints are excluded from affected models. Public case records contain only predefined categories, coded flags, and numeric elapsed times; case-level clinical prose is excluded. Symptom recovery, local lesion outcome, first post-flare response, and best index-course response are distinct outcomes.</p>
<p>Primary comparative analyses include site-reviewed suspected flare, including Possible cases. Higher-confidence and other prespecified sensitivities address classification and timing. Firth logistic models assess baseline and response associations. PFS and OS use time-dependent flare exposure, with day-30 and day-60 landmark sensitivities. Model populations, endpoint exclusions, covariate missingness, successful bootstrap counts, and qualifications accompany the estimates.</p>
<p>Published aggregate outputs regenerate the tables and figures. Generic analysis code is provided with synthetic inputs; confidential clinical source data are needed to rerun the empirical cohort analysis. Original clinical identifiers, clinical calendar dates, and case-level clinical prose are excluded. Synthetic inputs use artificial identifiers and dates.</p>
<p><a href='source.zip' download>Download source archive</a> <a href='Tables/Display_registry.csv' download>Download display registry</a></p>
""" + "".join(code_sections) + "</section>"
nav = "<a href='#pipeline'>Shared pipeline</a>" + "".join(f"<a href='#{name}'>{html.escape(name.replace('_', ' '))}</a>" for name in CAPTIONS) + "<a href='#references'>References</a>"
page = """<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width, initial-scale=1'>
<title>Tumor flare after bispecific antibody therapy supplementary analysis</title>
<style>
*{box-sizing:border-box}html{scroll-behavior:smooth}body{margin:0;color:#202830;background:#fff;font:16px/1.5 Arial,Helvetica,sans-serif}header{padding:38px max(24px,calc((100vw - 1280px)/2));background:#f0f3f5;border-bottom:1px solid #d9e0e4}h1{font-size:30px;line-height:1.2;max-width:1050px;margin:0 0 14px}header p{max-width:1050px;margin:0;color:#414e58}.layout{display:grid;grid-template-columns:190px minmax(0,1fr);gap:28px;max-width:1360px;margin:auto;padding:28px 24px}nav{position:sticky;top:20px;align-self:start;border-left:3px solid #d9e0e4;padding-left:14px}nav a{display:block;margin:8px 0;color:#1b5678;font-size:14px;text-decoration:none}main{min-width:0}section{margin-bottom:48px;padding-bottom:36px;border-bottom:1px solid #d9e0e4;scroll-margin-top:20px}h2{font-size:23px;line-height:1.25;margin:0 0 14px}p{margin:0 0 14px}.caption{font-size:15px;line-height:1.5}a{color:#1b5678;margin-right:12px}details{margin:14px 0;background:#f6f8f9;border:1px solid #d9e0e4;border-radius:4px}summary{cursor:pointer;padding:11px 14px;font-size:14px;font-weight:600}details p{padding:0 14px}pre{margin:0;padding:16px;overflow:auto;white-space:pre;tab-size:2;line-height:1.5;font:12px/1.5 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;background:#f7f8f9}figure{margin:18px 0}figure img{display:block;width:100%;height:auto;background:white}.table-scroll{overflow-x:auto;margin:18px 0 20px;border:1px solid #d9e0e4}table{border-collapse:collapse;width:100%;font-size:13px}th,td{border:1px solid #d9e0e4;padding:8px 10px;text-align:left;vertical-align:middle;white-space:normal;min-width:110px}th{background:#e9eff3;font-weight:600}tbody tr:nth-child(even){background:#f7f9fa}td:first-child,th:first-child{min-width:170px}footer{max-width:1360px;margin:auto;padding:0 24px 40px;font-size:13px;color:#4b5861}@media(max-width:850px){.layout{display:block;padding:20px 16px}nav{position:static;display:flex;flex-wrap:wrap;gap:4px 12px;margin-bottom:26px}nav a{display:inline-block}header{padding:24px 16px}h1{font-size:25px}h2{font-size:20px}}@media print{nav,details,a[download]{display:none}.layout{display:block;padding:0}header{background:white;padding:0 0 20px}section{break-inside:avoid}figure img{max-height:90vh;object-fit:contain}table{font-size:9px}}
</style></head><body><header><h1>Tumor flare after bispecific antibody therapy<br>Supplementary analysis and publication displays</h1><p>US Myeloma Immunotherapy Consortium. Tables and figures from the canonical analysis using site-reviewed EDC data as the authority for conflicting clinical findings.</p></header><div class='layout'><nav aria-label='Supplement contents'>""" + nav + "</nav><main>" + shared + "".join(sections) + reference_section + "</main></div><footer>Counts and estimates describe documented suspected flare in the available retrospective cohort. Clinical attribution, classification uncertainty, and missing data are reported with each analysis.</footer></body></html>"
(OUTPUT / "index.html").write_text(page, encoding="utf-8")
checks = {"tables": len(formatted_tables), "figures": len(FIGURE_FILES), "registry_rows": len(registry),
          "table1_groups": 2, "table2_groups": 3, "case_rows": len(formatted_tables["Table_S1"]),
          "html_self_contained_images": page.count("data:image/png;base64,"), "table_order": list(TABLE_BUILDERS),
          "figure_order": list(FIGURE_FILES), "generic_source_files": len(GENERIC_R_FILES), "rendered_count_cells_verified": rendered_count_cells,
          "doi_references": len(references), "html_sha256": sha256(OUTPUT / "index.html"), "caption_review": "Natural-writing review applied; uncertainty and attribution preserved"}
(OUTPUT / "publication_checks.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
print(json.dumps(checks, indent=2))
