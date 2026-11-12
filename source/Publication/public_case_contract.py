"""Closed contract for the public categorical case object.

Every field is a predefined category, a neutral label, or a numeric elapsed value.
Free text and unexpected fields are rejected before any display or archive build.
"""
import csv
import json
import math
import re
from pathlib import Path

NUMERIC_FIELDS = ("onset_days", "peak_days", "resolution_days", "onset_lower", "onset_upper")
BOOL_FIELDS = ("bridging", "comparative_included", "onset_estimated", "peak_estimated", "resolution_estimated",
               "first_response_on_index", "subsequent_response_on_index")
BINARY_FIELDS = ("pain_flare", "palpable_enlargement", "imaging_worsening", "new_lesion", "neurologic_cns",
                 "biochemical_or_organ", "paraskeletal_involvement", "cns_involvement", "spinal_involvement",
                 "organ_dysfunction", "flare_hospitalization", "icu_due_to_flare", "pathologic_fracture",
                 "cord_compression", "overlap_crs", "overlap_icans", "overlap_infection", "alternative_fracture",
                 "alternative_toxicity", "other_competing_explanation", "management_steroids",
                 "management_radiation", "management_surgery", "management_toci_anakinra", "true_pd_within_14d")
VOCABULARY = {
    "product": ("Teclistamab", "Talquetamab", "Elranatamab", "Linvoseltamab"),
    "confidence": ("Definite", "Probable", "Possible"),
    "symptom_status": ("RESOLVED", "IMPROVED_ONGOING", "ONGOING_NOT_IMPROVED", "DEATH_BEFORE_RESOLUTION", "UNKNOWN"),
    "local_outcome": ("RESOLVED", "IMPROVED", "STABLE", "WORSENED", "NOT_REASSESSED", "UNKNOWN"),
    "focality": ("SINGLE_SITE", "MULTIFOCAL", "DIFFUSE_SYSTEMIC", "UNKNOWN"),
    "known_or_new_site": ("KNOWN_SITE", "NEW_SITE", "UNCLEAR"),
    "skeletal_flare_type": ("LYTIC", "MIXED", "NONE", "OTHER_OSSEOUS", "UNKNOWN"),
    "true_emd_classification": ("NONE", "NONVISCERAL", "VISCERAL", "VISCERAL_AND_NONVISCERAL", "NOT_APPLICABLE", "UNKNOWN"),
    "bsab_action": ("CONTINUED", "DELAYED", "HELD_RESUMED", "DISCONTINUED", "UNKNOWN"),
    "first_postflare_response": ("sCR", "CR", "VGPR", "PR", "MR", "SD", "PD", "NE", "UNKNOWN", "NOT_APPLICABLE"),
    "subsequent_response": ("sCR", "CR", "VGPR", "PR", "MR", "SD", "PD", "NE", "NONE_DOCUMENTED", "NOT_REASSESSED", "UNKNOWN", "NOT_APPLICABLE"),
    "best_index_response": ("sCR", "CR", "VGPR", "PR", "MR", "SD", "PD", "UNAVAILABLE"),
    "alternative_diagnosis": ("NONE_IDENTIFIED", "TRUE_PROGRESSION", "UNCERTAIN"),
    "onset_timing_type": ("RECORDED_POINT", "ESTIMATED_POINT", "INTERVAL", "UNAVAILABLE"),
    "peak_timing_type": ("RECORDED_POINT", "ESTIMATED_POINT", "UNAVAILABLE"),
    "resolution_timing_type": ("RECORDED_POINT", "ESTIMATED_POINT", "UNAVAILABLE"),
    "first_response_category_estimated": ("YES", "NO", "NOT_REPORTED"),
    "later_local_assessment_context": ("AFTER_RADIATION", "AFTER_SALVAGE", "AFTER_CART", "AFTER_RADIATION_AND_SALVAGE", "NOT_REPORTED"),
    "later_local_after_radiation": ("YES", "NO", "NOT_REPORTED"),
    "later_local_after_salvage": ("YES", "NO", "NOT_REPORTED"),
    "later_local_after_cart": ("YES", "NO", "NOT_REPORTED"),
    "later_local_separate_from_systemic_response": ("YES", "NOT_REPORTED"),
    "peak_before_point_onset": ("YES", "NO", "NOT_EVALUABLE"),
}
VOCABULARY.update({field: ("TRUE", "FALSE", "UNKNOWN") for field in BOOL_FIELDS})
VOCABULARY.update({field: ("YES", "NO", "UNKNOWN", "NOT_APPLICABLE") for field in BINARY_FIELDS})
FIELD_ORDER = ("case_label", "product", "confidence", "bridging", "comparative_included", *NUMERIC_FIELDS,
               "onset_estimated", "peak_estimated", "resolution_estimated", "onset_timing_type", "peak_timing_type",
               "resolution_timing_type", "symptom_status", "local_outcome", *BINARY_FIELDS,
               "focality", "known_or_new_site", "skeletal_flare_type", "true_emd_classification", "bsab_action",
               "first_postflare_response", "first_response_on_index", "first_response_category_estimated",
               "subsequent_response", "subsequent_response_on_index", "best_index_response", "alternative_diagnosis",
               "later_local_assessment_context", "later_local_after_radiation", "later_local_after_salvage",
               "later_local_after_cart", "later_local_separate_from_systemic_response", "peak_before_point_onset")
if len(FIELD_ORDER) != len(set(FIELD_ORDER)):
    raise RuntimeError("Duplicate public case contract fields")

def validate_cases(rows):
    if not rows:
        raise ValueError("Public case object is empty")
    labels = set()
    for row in rows:
        if set(row) != set(FIELD_ORDER):
            raise ValueError(f"Public case columns differ from closed contract: {set(row) ^ set(FIELD_ORDER)}")
        label = row["case_label"]
        if not re.fullmatch(r"TF[0-9]{2}", label) or label in labels:
            raise ValueError("Neutral public case labels must be unique TFnn values")
        labels.add(label)
        for field, allowed in VOCABULARY.items():
            if row[field] not in allowed:
                raise ValueError(f"Public case value outside finite vocabulary: {field}")
        for field in NUMERIC_FIELDS:
            value = row[field]
            if value == "":
                continue
            if not re.fullmatch(r"[0-9]+(?:\.[0-9]+)?", value) or not math.isfinite(float(value)) or float(value) > 36600:
                raise ValueError(f"Public chronology must be a finite nonnegative elapsed value: {field}")
        for event in ("onset", "peak", "resolution"):
            point = row[event + "_days"]
            timing = row[event + "_timing_type"]
            if timing in ("RECORDED_POINT", "ESTIMATED_POINT") and not point:
                raise ValueError("Point timing category requires an elapsed value")
            if timing in ("UNAVAILABLE", "INTERVAL") and point:
                raise ValueError("Unknown or interval timing cannot contain an inferred point")
            if timing == "ESTIMATED_POINT" and row[event + "_estimated"] != "TRUE":
                raise ValueError("Estimated timing must carry the corresponding flag")
            if timing == "RECORDED_POINT" and row[event + "_estimated"] != "FALSE":
                raise ValueError("Recorded timing must carry the corresponding flag")
        if row["onset_timing_type"] == "INTERVAL" and (not row["onset_lower"] or not row["onset_upper"] or float(row["onset_lower"]) > float(row["onset_upper"])):
            raise ValueError("Onset interval requires ordered elapsed bounds")
        if row["onset_timing_type"] != "INTERVAL" and (row["onset_lower"] or row["onset_upper"]):
            raise ValueError("Onset bounds require an interval timing category")
        contexts = {
            "NOT_REPORTED": ("NOT_REPORTED", "NOT_REPORTED", "NOT_REPORTED", "NOT_REPORTED"),
            "AFTER_RADIATION": ("YES", "NO", "NO", "YES"),
            "AFTER_SALVAGE": ("NO", "YES", "NO", "YES"),
            "AFTER_CART": ("NO", "NO", "YES", "YES"),
            "AFTER_RADIATION_AND_SALVAGE": ("YES", "YES", "NO", "YES"),
        }
        if tuple(row[field] for field in ("later_local_after_radiation", "later_local_after_salvage", "later_local_after_cart", "later_local_separate_from_systemic_response")) != contexts[row["later_local_assessment_context"]]:
            raise ValueError("Local-assessment treatment flags contradict the coded context")
        expected_peak = ("YES" if float(row["peak_days"]) < float(row["onset_days"]) else "NO") if row["peak_days"] and row["onset_days"] else "NOT_EVALUABLE"
        if row["peak_before_point_onset"] != expected_peak:
            raise ValueError("Chronology flag contradicts elapsed values")
    return rows

def read_coded_cases(path):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        return validate_cases(list(csv.DictReader(handle)))

def public_generic_source(filename, source):
    """Suppress narrative exports in the public copy without changing models."""
    if filename != "summary_functions.R":
        return source
    pattern = r'  for\(field in intersect\(c\("first_response_on_index".*?  tf_write\(summary,out,"case_summary"\)'
    replacement = ('  for(field in intersect(c("first_response_on_index","subsequent_response_on_index"),names(c)))summary[[field]]<-c[[field]]\n'
                   '  tf_write(summary,out,"case_summary")')
    result, count = re.subn(pattern, replacement, source, flags=re.S)
    if count != 1:
        raise ValueError("Public summary export suppression did not match its expected code boundary")
    return result

def write_contract(path):
    contract = {"schema": "tumor-flare-public-case-coded-v1", "fields": list(FIELD_ORDER),
                "neutral_label_pattern": "TF[0-9]{2}", "numeric_elapsed_fields": list(NUMERIC_FIELDS),
                "numeric_unit": "days after first step-up dose", "empty_numeric_value": "unavailable",
                "finite_category_vocabulary": {key: list(value) for key, value in VOCABULARY.items()},
                "local_assessment_not_reported": "Absence of a reported additional-treatment assessment does not establish absence of that treatment.",
                "response_attribution": "First and subsequent on-index flags qualify the corresponding response categories. Additional-treatment local assessments are separate from index-course systemic response.",
                "estimated_response_flag": "YES identifies a source-supported viewer-estimated response category. NOT_REPORTED means that estimation status is not established."}
    Path(path).write_text(json.dumps(contract, indent=2), encoding="utf-8")
