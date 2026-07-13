# Changelog

All notable changes to the analysis code and aggregate outputs in this
repository are documented here. The dissertation document is a separate
deliverable and is not tracked in this repository.

## v1.0.3 — BFDP credibility analysis, Figure 4.3, and Appendix A figures

- Added `scripts/bfdp_analysis.R`: Wakefield (2007) Bayesian false-discovery
  probability for every DCIS gene association, complementing the Benjamini–Hochberg
  FDR, with a 3/5/8 prior-odds-ratio sensitivity sweep. Five associations are
  noteworthy under the primary prior (*ATM* truncating BFDP < 0.001, *CHEK2*
  truncating 0.010, *CHEK2* missense 0.028, *BRCA2* truncating 0.035, *TP53*
  missense 0.046); *TP53* missense is prior-sensitive (rises to 0.071 under the
  tightest prior). Outputs `outputs/tables/BFDP_analysis_DCIS.csv`,
  `BFDP_prior_sensitivity.csv`, `Table_4_3a_BFDP.csv`.
- Added `scripts/figure_4_3_volcano.R`: BFDP-encoded volcano of the DCIS
  associations (x = log2 OR, y = −log10 p, colour = BFDP, shape = variant class)
  → `outputs/figures/figure_4_3_volcano.png` (+ vector PDF).
- Added `scripts/appendix_figures_A1_A15.R`: Appendix A supplementary figures
  A1–A15 (per-outcome and per-subtype forest plots, cohort QC, ascertainment
  sensitivity, leave-one-study-out, carrier depletion) → `outputs/figures/appendix/`.
- Figure hygiene: long panel titles wrapped, *A11* facet labels no longer clipped,
  *A2* uses the Unicode ≥ sign; all figures re-rendered as PNG + vector PDF.
- Positive-control specification note (aggregate outputs only): the invasive
  ever-HRT estimate is reported under the fuller both-arm, study-adjusted
  specification (OR 1.09, increased risk, as expected) rather than the simpler
  unadjusted model (OR 0.88); post-menopausal status is reported relative to
  pre-menopausal women. Neither is used as a validation benchmark.

## v1.0.2 — KARMA-updated BBD sensitivity

- Added `scripts/karma_bbd_sensitivity.R`: robustness of the BBD truncating-variant
  results to the KARMA study's updated BBD-status data, comparing the primary,
  KARMA-excluded (drop-all) and KARMA-updated cohorts across the nine genes and the
  overall / non-proliferative / proliferative panels. *BRCA2* overall remains the
  only nominally significant, stable association (OR 2.87); no gene reaches FDR
  significance in any panel or version.
- The KARMA update and individual-level data are not distributed (data agreement);
  only the code is released. Aggregate cells with <5 carriers are non-disclosive
  and are not included.

## v1.0.1 — License and DOI

- Added the MIT License (© 2026 Zixi Yao).
- Recorded the Zenodo concept DOI (all versions):
  [10.5281/zenodo.21280365](https://doi.org/10.5281/zenodo.21280365).

## v1.0.0 — Initial public release

- First public, citable release of the analysis code and non-disclosive
  aggregate outputs for the dissertation *Genetic Associations Between Benign
  Breast Disease, Carcinoma in Situ and Breast Cancer* (MPhil Population Health
  Sciences, University of Cambridge, 2025–26).
- Analysis: Firth penalised-likelihood models for rare truncating and missense
  variants across nine breast-cancer susceptibility genes and three outcomes
  (BBD, DCIS, LCIS); hormonal and reproductive factor models; ascertainment and
  confidence-interval-method sensitivity analyses.
- Individual-level BRIDGES/BCAC data and the BCAC data dictionary are excluded
  under the data-access agreement (see README, "Data and code availability").
- Absolute local paths removed from scripts; scripts are run from the
  repository root.
