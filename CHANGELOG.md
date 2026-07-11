# Changelog

All notable changes to the analysis code and aggregate outputs in this
repository are documented here. The dissertation document is a separate
deliverable and is not tracked in this repository.

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
