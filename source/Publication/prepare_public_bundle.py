#!/usr/bin/env python3
"""Stage an explicit public allowlist. Does not deploy or publish.

Usage: python prepare_public_bundle.py ANALYSIS_DIR PUBLICATION_DIR CODE_DIR BUNDLE_DIR
"""
from __future__ import annotations
import csv
import hashlib
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path
from public_case_contract import read_coded_cases, public_generic_source

analysis, publication, code, bundle = map(Path, sys.argv[1:5])
here = Path(__file__).resolve().parent
bundle.mkdir(parents=True, exist_ok=True)
if (bundle / ".git").exists():
    raise RuntimeError("Stage the public bundle outside an existing Git repository")
if (analysis / "case_summary.csv").exists():
    raise RuntimeError("Public staging requires coded inputs; clinical case_summary.csv is forbidden")
read_coded_cases(analysis / "case_coded.csv")

forbidden_names = re.compile(r"(^|_)(course_id|patient_id|mrn|source_reference|narrative|calendar_date|start_date|description|evidence|date_notes|source_path)($|_)", re.I)
date_pattern = re.compile(r"\b(?:19|20)\d{2}[-/]\d{1,2}[-/]\d{1,2}\b")
clinical_id_pattern = re.compile(r"\b(?:CCF|TEC|TALQ|MDACC)[_-]?\d{2,}\b", re.I)
allowed_paths: set[str] = set()

def validate_csv(path: Path, synthetic: bool = False) -> None:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        bad = [name for name in (reader.fieldnames or []) if forbidden_names.search(name)]
        if bad and not synthetic:
            raise ValueError(f"Private source columns in {path.name}: {bad}")
        for row in reader:
            for key, value in row.items():
                value = value or ""
                if date_pattern.search(value) and key not in ("run_date",) and not synthetic:
                    raise ValueError(f"Calendar date in {path.name}:{key}")
                if clinical_id_pattern.search(value):
                    raise ValueError(f"Original clinical identifier in {path.name}:{key}")
                if "/Users/" in value or "/private/" in value:
                    raise ValueError(f"Private absolute path in {path.name}:{key}")

def copy(source: Path, relative: str, synthetic: bool = False) -> None:
    if not source.exists():
        raise FileNotFoundError(source)
    if source.suffix.lower() == ".csv":
        validate_csv(source, synthetic=synthetic)
        if source.name in ("Table_S1.csv", "case_coded.csv"):
            read_coded_cases(source)
    destination = bundle / relative
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, destination)
    allowed_paths.add(relative)

figure_stems = (
    "Figure_1_Cohort_flow", "Figure_2_Onset_and_chronology", "Figure_3_Clinical_spectrum",
    "Figure_4_Association_models", "Figure_5_Recognition_framework", "Figure_S1_Landmark_survival")
for stem in figure_stems:
    for suffix in (".png", ".pdf"):
        copy(publication / "Figures" / (stem + suffix), "assets/" + stem + suffix)
for stem in ("Table_1", "Table_2", "Table_S1", "Table_S2", "Table_S3", "Table_S4", "Table_S5", "Display_registry"):
    copy(publication / "Tables" / (stem + ".csv"), "tables/" + stem + ".csv")
registry_path = bundle / "tables/Display_registry.csv"
with registry_path.open(newline="", encoding="utf-8") as handle:
    registry_reader = csv.DictReader(handle)
    registry_fields = registry_reader.fieldnames
    registry_rows = list(registry_reader)
for row in registry_rows:
    row["dataset"] = row["dataset"].replace("Tables/", "tables/").replace("Display_data/", "tables/")
with registry_path.open("w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, registry_fields)
    writer.writeheader(); writer.writerows(registry_rows)
for stem in ("Figure_1", "Figure_3", "Figure_4", "Figure_5", "Figure_S1_KM", "Figure_S1_Risk"):
    copy(publication / "Display_data" / (stem + ".csv"), "tables/" + stem + ".csv")

index = (publication / "index.html").read_text(encoding="utf-8")
index = index.replace("href='Tables/", "href='tables/").replace("href='Figures/", "href='assets/")
(bundle / "index.html").write_text(index, encoding="utf-8")
(bundle / ".nojekyll").write_text("", encoding="utf-8")
allowed_paths.update(("index.html", ".nojekyll"))

source_dir = bundle / "source"
source_dir.mkdir(exist_ok=True)
GENERIC_R_FILES = ("run_analysis.R", "analysis_functions.R", "survival_functions.R",
                   "reconciliation_functions.R", "model_functions.R", "summary_functions.R", "ground_truth_tests.R")
code_files = [code / name for name in GENERIC_R_FILES]
for path in code_files:
    source_text = path.read_text(encoding="utf-8")
    if "/Users/" in source_text or clinical_id_pattern.search(source_text):
        raise ValueError(f"Private identifier or path in code: {path.name}")
    public_source = public_generic_source(path.name, source_text)
    destination = bundle / "source/Code" / path.name
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(public_source, encoding="utf-8")
    allowed_paths.add("source/Code/" + path.name)
copy(code / "model_seed_manifest.csv", "source/Code/model_seed_manifest.csv")
for filename in ("build_figures.R", "build_publication.py", "prepare_public_bundle.py", "public_case_contract.py"):
    copy(here / filename, "source/Publication/" + filename)

# These files contain reported aggregate results or the curated, neutral Table S1
# fields already present in the public manuscript supplement. Raw case timelines,
# survival observations, endpoint exclusion rows, and original source logs are absent.
reported_outputs = (
    "table1_baseline", "table2_clinical_course", "flow", "coverage", "incidence",
    "presentation", "anatomy_binary", "timing_summary", "case_coded", "models",
    "survival_models", "clinical_comparisons", "response_distribution", "response_thresholds",
    "hospitalization", "missingness", "bootstrap_diagnostics", "landmark_models",
    "landmark_km_steps", "landmark_risk_counts")
for name in reported_outputs:
    copy(analysis / (name + ".csv"), "source/data/reported/" + name + ".csv")
copy(analysis / "public_case_vocabulary.json", "source/data/reported/public_case_vocabulary.json")

dictionary = code.parent / "Data/controlled_input_dictionary.csv"
if not dictionary.exists():
    raise FileNotFoundError("Controlled input dictionary required")
copy(dictionary, "source/Config/controlled_input_dictionary.csv")
demo_dir = code.parent / "Data/Demo"
if not demo_dir.exists():
    raise FileNotFoundError("Synthetic demonstration inputs required")
synthetic_marker = (demo_dir / "SYNTHETIC.txt").read_text(encoding="utf-8")
if "fully synthetic" not in synthetic_marker.lower():
    raise ValueError("Synthetic input attestation required")
for filename in ("master_courses.csv", "edc_cases.csv"):
    with (demo_dir / filename).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if row.get("course_id") and not re.fullmatch(r"D\d{4}", row["course_id"]):
                raise ValueError("Synthetic course ID namespace violated")
            if row.get("patient_cluster_id") and not re.fullmatch(r"DP\d{4}", row["patient_cluster_id"]):
                raise ValueError("Synthetic cluster ID namespace violated")
            if row.get("site_key") and not row["site_key"].startswith("demo_"):
                raise ValueError("Synthetic site namespace violated")
for filename in ("master_courses.csv", "edc_cases.csv", "verified_corrections.csv"):
    copy(demo_dir / filename, "source/Demo/" + filename, synthetic=True)
copy(demo_dir / "SYNTHETIC.txt", "source/Demo/SYNTHETIC.txt")
copy(publication / "reference_metadata_public.json", "source/References/reference_metadata_public.json")
copy(here / "requirements-publication.txt", "source/requirements-publication.txt")
copy(here / "publication_runtime.json", "source/publication_runtime.json")
lock = code.parent / "renv.lock"
if not lock.exists():
    lock = code / "renv.lock"
if not lock.exists():
    raise FileNotFoundError("The reproducibility bundle requires renv.lock")
copy(lock, "source/renv.lock")

readme = """# Tumor flare manuscript reproducibility source

This archive contains the generic cohort-analysis code, publication code, reported
aggregate outputs, and synthetic demonstration inputs for the US Myeloma
Immunotherapy Consortium tumor-flare study.

Site-reviewed EDC findings override conflicting master-dataset findings for the
same clinical concept. Unknown authoritative findings remain unknown. The code
keeps first post-flare response, best index-course response, symptom recovery, and
local lesion outcome separate and uses recorded flare onset in survival models.

## Rebuild publication displays

Publication displays consume `data/reported/case_coded.csv`, its closed category
vocabulary, and the reported aggregate CSV outputs. Clinical `case_summary.csv`
and its qualifications are not publication inputs.

From this folder, after installing the R dependencies recorded in `renv.lock`:

```sh
python3 -m pip install -r requirements-publication.txt
Rscript Publication/build_figures.R data/reported publication
python3 Publication/build_publication.py data/reported publication Code
```

The Python publication renderer requires Pillow for embedded figure previews.
R figures use Helvetica in vector PDF and 320 dpi PNG exports. A conceptual
framework can also be regenerated without clinical data:

```sh
Rscript Publication/build_figures.R data/reported publication --framework-only
```

## Cohort analysis

The generic runner and controlled input dictionary describe the input contract.
Run the deterministic synthetic demonstration from this folder:

```sh
Rscript Code/run_analysis.R Demo demo_results 20
```

This small draw count exercises data reconciliation and the model workflow;
it is not used for empirical manuscript inference. The full manuscript analysis
uses 1,000 patient-cluster bootstrap draws per logistic model. The actual master
workbook, completed EDCs, original patient/course identifiers, calendar dates,
clinical narratives, and patient-level survival observations are not distributed.
Therefore the public archive can reproduce publication displays from reported
outputs and demonstrate the analysis pipeline; recomputing the empirical cohort
estimates requires authorized access to the original clinical data.

The public case table contains only neutral labels, predefined coded categories
and flags, and numeric elapsed times. Every value is checked against a closed
contract. All case-level qualification, chronology-review, and local-assessment
prose is excluded. The full clinical qualifications remain in the local,
controlled case_summary.csv and are not part of this archive. Raw case timelines
and reconciliation histories are excluded. Kaplan-Meier data contain aggregate survival steps and
risk counts, not individual follow-up records.

The public copy of summary_functions.R suppresses all case narrative export
fields and chronology prose templates. This changes descriptive case export
columns only; the reconciliation, response derivation, and statistical models
are unchanged. The local clinical source modules are not modified.

## Checks and provenance

`tables/Display_registry.csv` in the site package maps every display to its
analysis outputs, executable publication module, source hashes, and checks.
`source_manifest.json` records SHA-256 hashes for this archive's files. The canonical
analysis exports report model populations, confidence intervals, successful
bootstrap counts, and estimation qualifications. No clear association is not
interpreted as evidence of equivalent outcomes.
"""
(source_dir / "README.md").write_text(readme, encoding="utf-8")
allowed_paths.add("source/README.md")

manifest = []
for path in sorted(source_dir.rglob("*")):
    if path.is_file() and path.name != "source_manifest.json":
        manifest.append({"path": path.relative_to(source_dir).as_posix(), "sha256": hashlib.sha256(path.read_bytes()).hexdigest(), "bytes": path.stat().st_size})
(source_dir / "source_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
allowed_paths.add("source/source_manifest.json")
with zipfile.ZipFile(bundle / "source.zip", "w", compression=zipfile.ZIP_DEFLATED) as archive:
    for path in sorted(source_dir.rglob("*")):
        if path.is_file():
            archive.write(path, "source/" + path.relative_to(source_dir).as_posix())
allowed_paths.add("source.zip")
unexpected = [p.relative_to(bundle).as_posix() for p in bundle.rglob("*") if p.is_file() and p.relative_to(bundle).as_posix() not in allowed_paths]
if unexpected:
    raise ValueError(f"Unallowlisted files in bundle: {unexpected}")
shutil.copyfile(bundle / "source.zip", publication / "source.zip")
checks = {"staged_only": True, "figure_files": len(figure_stems)*2, "table_files": 8,
          "source_files": len(manifest)+1, "reported_output_files": len(reported_outputs),
          "raw_clinical_inputs_included": False, "curated_neutral_case_table_included": True,
          "case_values_closed_vocabulary": True, "case_free_prose_included": False,
          "raw_landmark_rows_included": False, "raw_case_timeline_included": False,
          "source_zip_sha256": hashlib.sha256((bundle / "source.zip").read_bytes()).hexdigest()}
(publication / "public_bundle_checks.json").write_text(json.dumps(checks, indent=2), encoding="utf-8")
print(json.dumps(checks, indent=2))
