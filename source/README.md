# Tumor flare manuscript reproducibility source

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
