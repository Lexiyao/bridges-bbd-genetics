# Genetic Associations Between Benign Breast Disease, Carcinoma in Situ and Breast Cancer

MPhil in Population Health Sciences, University of Cambridge, 2025–26

## Project Overview

This project examines associations between pathogenic variants in nine breast cancer susceptibility genes (BRCA1, BRCA2, PALB2, CHEK2, ATM, BARD1, RAD51C, RAD51D, TP53) and three outcomes: benign breast disease (BBD), ductal carcinoma in situ (DCIS), and lobular carcinoma in situ (LCIS). Hormonal and reproductive risk factors are examined as secondary objectives.

Data source: BRIDGES/BCAC consortium, phenotype file v17.

## Citable archive

This repository contains the analysis code and aggregate (non-disclosive)
outputs for the dissertation. A citable snapshot is archived on Zenodo:

> DOI: *to be assigned on release* — cite the "all versions" (concept) DOI.

The dissertation document itself is a separate deliverable and is not
redistributed here. See the [CHANGELOG](CHANGELOG.md) for version history.

## Directory Structure

```
project/
├── scripts/                        # All R analysis scripts
│   ├── BBD_truncating_analysis_v3.R       # FINAL: BBD truncating variant analysis (Firth)
│   ├── DCIS_LCIS_truncating_analysis.R    # FINAL: DCIS/LCIS truncating variant analysis (Firth)
│   ├── missense_analysis_corrected.R      # FINAL: Missense variant analysis (Firth, corrected)
│   ├── objective2_hormonal_analysis.R     # FINAL: Hormonal/reproductive factor analysis (glm)
│   ├── export_all_results.R               # Exports combined results to Excel
│   ├── BBD_firth_patch.R                  # Firth regression patch/helper
│   └── 第一个部分.R                        # Early exploratory script (archived)
├── outputs/
│   ├── figures/                    # Forest plots (PNG, 300 dpi)
│   │   ├── forest_BBD_overall.png
│   │   ├── forest_BBD_nonproliferative.png
│   │   ├── forest_BBD_proliferative.png
│   │   ├── forest_DCIS_truncating.png
│   │   ├── forest_LCIS_truncating.png
│   │   ├── forest_missense_BBD_corrected.png
│   │   ├── forest_missense_BBD_nonprol_corrected.png
│   │   ├── forest_missense_BBD_prol_corrected.png
│   │   ├── forest_missense_DCIS_corrected.png
│   │   ├── forest_missense_LCIS_corrected.png
│   │   ├── forest_hormonal_corrected.png
│   │   ├── forest_DCIS_significant_summary.png   # Headline: 5 FDR-sig DCIS associations
│   │   ├── heatmap_overview.png                  # 9 genes × 6 outcome/variant cells
│   │   ├── gradient_stage_truncating.png         # OR by stage (BBD → DCIS) for ATM/CHEK2/BRCA2
│   │   └── participant_flow.png                  # CONSORT-style cohort derivation
│   ├── tables/                     # Result Excel files (FINAL versions)
│   │   ├── BBD_truncating_FINAL.xlsx
│   │   ├── DCIS_LCIS_truncating_results.xlsx
│   │   ├── missense_results_FINAL.xlsx
│   │   └── objective2_hormonal_corrected.xlsx
│   └── session_info/               # R session info for reproducibility
├── CHANGELOG.md
└── README.md
```

**Raw BRIDGES data files are NOT stored here** (data governance requirement).
Place them in the project root before running scripts:
- `concept_807_zhang_bridges_pheno_v17.txt`
- `concept_807_zhang_bridges_truncating.csv`
- `concept_807_zhang_bridges_missense.csv`

## Analysis Pipeline

```
raw data → data_QC.R → analysis scripts → outputs/tables/ + outputs/figures/
             ↓
         outputs/QC/   (QC report + figures, review before proceeding)
```

Run in order:
1. `scripts/data_QC.R` — **Run first**: data integrity, phenotype QC, carrier frequency checks → `outputs/QC/`
2. `scripts/BBD_truncating_analysis_v3.R` — BBD overall and subtype truncating models → `BBD_truncating_FINAL.xlsx`
3. `scripts/DCIS_LCIS_truncating_analysis.R` — DCIS/LCIS truncating models → `DCIS_LCIS_truncating_results.xlsx`
4. `scripts/missense_analysis_corrected.R` — All missense models (BBD, DCIS, LCIS) → `missense_results_FINAL.xlsx`
5. `scripts/objective2_hormonal_analysis.R` — Hormonal/reproductive factors → `objective2_hormonal_corrected.xlsx`
6. `scripts/sensitivity_familial_population.R` — ascertainment sensitivity (all-studies vs population-only vs famHist-adjusted) → `sensitivity_familial_population.xlsx`
7. `scripts/export_all_results.R` — (optional) combine results if running from a joint session

**Verification / audit (independent of the analysis, run any time):**
- `scripts/verify_carrier_counts.R` — re-counts every Carriers(cases)/Carriers(controls) value straight from the raw genotype files and compares to the result tables → `outputs/QC/carrier_count_verification.xlsx` (72/72 confirmed, 0 mismatch).
- `scripts/compare_ci_methods.R` — confirms Wald vs profile-likelihood CIs are materially identical for the significant DCIS genes → `outputs/tables/ci_method_check.xlsx`.

## Computational Environment

- R 4.5.2 (2025-10-31), macOS 14 (Sonoma)
- Key packages: tidyverse, logistf, writexl, ggplot2, data.table
- Full package versions: see `outputs/session_info/` after running scripts

## Data and code availability

The analysis code in this repository is openly available and archived with a
persistent identifier on Zenodo (DOI: *to be assigned on release*).

Individual-level BRIDGES phenotype and genotype data were accessed under an
approved BCAC secondary-use agreement and **cannot be redistributed**. Raw
participant-level data are not committed to version control and are excluded by
`.gitignore`. Bona fide researchers may apply for data access through the Breast
Cancer Association Consortium (http://bcac.ccge.medschl.cam.ac.uk/).

The BRIDGES data dictionary is BCAC-controlled documentation and is likewise not
redistributed here; approved users obtain it from BCAC.

All outputs committed here are non-disclosive summary statistics (aggregate
carrier counts, odds ratios and figures).

## Key Results (verified)

### Objective 1 — Genetic associations

| Outcome | Gene | Variant class | OR (95% CI) | FDR q |
|---------|------|--------------|-------------|-------|
| DCIS | ATM | Truncating | 3.99 (2.13–7.02) | <0.001 |
| DCIS | CHEK2 | Truncating | 2.66 (1.51–4.42) | 0.005 |
| DCIS | BRCA2 | Truncating | 2.82 (1.43–5.17) | 0.012 |
| DCIS | CHEK2 | Missense (CADD phred ≥20) | 2.01 (1.29–3.00) | 0.023 |
| DCIS | TP53 | Missense (CADD phred ≥20) | 3.65 (1.51–7.76) | 0.027 |
| BBD | — | — | No FDR-significant associations | — |
| LCIS | — | — | No FDR-significant associations | — |

### Objective 2 — Hormonal and reproductive factors

| Exposure | Outcome | OR per unit | FDR q |
|----------|---------|------------|-------|
| Parity (per birth) | BBD | 0.938 | 0.009 |
| Parity (per birth) | DCIS | 0.909 | 0.002 |
| Parity (per birth) | LCIS | 0.754 | 0.009 |
| Post-menopausal status | DCIS | 0.513 | <0.001 |
| Ever HRT use | BBD | 1.636 | <0.001 |

Sample sizes: BBD 1,765 cases / 12,715 controls; DCIS 1,663 / 34,424; LCIS 139 / 34,424.
All models: Firth's logistic regression (Obj 1) or standard glm (Obj 2), adjusted for age and study.
FDR: Benjamini-Hochberg across 9 genes (Obj 1) or 15 exposure-outcome pairs (Obj 2).

## Key Analysis Decisions

| Decision | Rationale |
|----------|-----------|
| Firth's penalised logistic regression | Handles complete/quasi-complete separation from rare variant carriers |
| FDR: Benjamini-Hochberg across 9 genes | Balances Type I error control with power in multi-gene exploration |
| European ancestry restriction (ethnicityClass == 1) | Reduces population stratification; PCs not available in released phenotype file |
| Complete case analysis for hormonal variables | MNAR by study centre; MI assumptions not met |
| No PC adjustment | PC columns absent from BRIDGES phenotype file v17; limitation documented in Section 5.4 |

## License

The code in this repository is released under the [MIT License](LICENSE)
(© 2026 Zixi Yao). Individual-level BRIDGES/BCAC data are not covered by this
license and remain controlled-access (see "Data and code availability").
