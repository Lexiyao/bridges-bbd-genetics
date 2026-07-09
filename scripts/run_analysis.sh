#!/usr/bin/env bash
# =============================================================================
# Reproduce the analysis and figures for
#   "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ
#    and Breast Cancer" — MPhil Population Health Sciences, Cambridge 2025-26
#
# Requires the raw BRIDGES data files in the repository root (NOT distributed
# here — see README, "Data and code availability"):
#   concept_807_zhang_bridges_pheno_v17.txt
#   concept_807_zhang_bridges_truncating.csv
#   concept_807_zhang_bridges_missense.csv
#
# Usage:  ./scripts/run_analysis.sh
# =============================================================================
set -euo pipefail
cd "$(dirname "$0")/.."                     # repository root, regardless of caller

echo "==> [1] data QC (run first: integrity, phenotype QC, carrier checks)"
Rscript scripts/data_QC.R

echo "==> [2] genetic and hormonal models (Firth / glm)"
for s in DCIS_LCIS_truncating_analysis BBD_truncating_analysis_v3 \
         missense_analysis_corrected objective2_hormonal_analysis \
         er_descriptive sensitivity_familial_population; do
  echo "    Rscript scripts/$s.R"
  Rscript "scripts/$s.R"
done

echo "==> [3] figures from the verified result tables"
Rscript scripts/build_exhibits.R    # forest plots, overview heatmap, headline summary
Rscript scripts/build_flowchart.R   # participant-flow diagram
Rscript scripts/build_dag.R         # collider-bias schematic
Rscript scripts/build_gradient.R    # OR-by-stage gradient (ATM/CHEK2/BRCA2)

echo "==> done. See outputs/tables/ and outputs/figures/."
