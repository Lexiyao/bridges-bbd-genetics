# =============================================================================
# VERIFICATION SCRIPT — independently reproduces all key thesis numbers
# =============================================================================

suppressMessages({ library(data.table); library(tidyverse); library(logistf) })

PASS <- 0L; FAIL <- 0L
chk <- function(label, got, expected, tol=0) {
  if (length(got)==0 || is.na(got)) {
    cat(sprintf("  FAIL <<<<< %-52s got=NA  expected=%s\n", label, expected))
    FAIL <<- FAIL+1L; return(invisible(NULL))
  }
  ok <- abs(as.numeric(got) - as.numeric(expected)) <= tol
  cat(sprintf("%s  %-52s got=%-10s expected=%s\n",
              if(ok) "  PASS" else "  FAIL <<<<<", label, round(as.numeric(got),3), expected))
  if (ok) PASS <<- PASS+1L else FAIL <<- FAIL+1L
  invisible(NULL)
}

fit_or <- function(data, gene_col) {
  cat(sprintf("    fitting %s ...", gene_col))
  md <- data |>
    select(outcome, carrier = all_of(gene_col), ageInt, study) |>
    mutate(study = factor(study)) |>
    drop_na()
  cat(sprintf(" n=%d, carriers=%d\n", nrow(md), sum(md$carrier)))
  fit <- logistf(outcome ~ carrier + ageInt + study,
                 data = md, firth = TRUE, pl = FALSE)   # pl=FALSE for speed
  i <- which(names(coef(fit)) == "carrier")
  # Wald CI when pl=FALSE
  se <- sqrt(diag(vcov(fit)))[i]
  list(OR  = exp(coef(fit)[i]),
       lo  = exp(coef(fit)[i] - 1.96*se),
       hi  = exp(coef(fit)[i] + 1.96*se),
       p   = fit$prob[i])
}

# ── 1. LOAD ────────────────────────────────────────────────────────────────────
cat("── Loading data ──\n")
ph <- fread("concept_807_zhang_bridges_pheno_v17.txt",
            na.strings=c("","NA","888","777","999")) |> as_tibble()
tr <- fread("concept_807_zhang_bridges_truncating.csv",
            na.strings=c("","NA")) |> as_tibble() |> select(-any_of(c("...1","V1")))
ms <- fread("concept_807_zhang_bridges_missense.csv",
            na.strings=c("","NA")) |> as_tibble() |> select(-any_of(c("...1","V1")))

GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")

# ── 2. BUILD CARRIER COLUMNS ──────────────────────────────────────────────────
tr_carrier <- tr |> select(BRIDGES_ID)
ms_carrier <- ms |> select(BRIDGES_ID)
for (g in GENES) {
  tc <- paste0(g, "_truncating");     if (tc %in% names(tr)) tr_carrier[[g]]         <- replace_na(as.integer(tr[[tc]]),  0L)
  mc <- paste0(g, "_CADD.phred.01"); if (mc %in% names(ms)) ms_carrier[[paste0(g,"_m")]] <- replace_na(as.integer(ms[[mc]]), 0L)
}
cat("Truncating cols built:", paste(intersect(GENES, names(tr_carrier)), collapse=", "), "\n")
cat("Missense cols built:  ", paste(paste0(GENES,"_m")[paste0(GENES,"_m") %in% names(ms_carrier)], collapse=", "), "\n\n")

dat <- ph |>
  left_join(tr_carrier, by="BRIDGES_ID") |>
  left_join(ms_carrier, by="BRIDGES_ID")

# ── 3. COHORTS ────────────────────────────────────────────────────────────────
eur0     <- dat |> filter(ethnicityClass==1, status==0, !is.na(ageInt))
bbd_all  <- eur0 |> filter(!is.na(BBD_history)) |> mutate(outcome=as.integer(BBD_history==1))
bbd_np   <- bbd_all |> filter(outcome==0 | (!is.na(BBD_type1) & BBD_type1==1)) |>
            mutate(outcome=if_else(!is.na(BBD_type1) & BBD_type1==1, 1L, 0L))
bbd_pr   <- bbd_all |> filter(outcome==0 | (!is.na(BBD_type1) & BBD_type1==2)) |>
            mutate(outcome=if_else(!is.na(BBD_type1) & BBD_type1==2, 1L, 0L))
ctrl     <- eur0 |> mutate(outcome=0L)
dcis_c   <- dat |> filter(ethnicityClass==1, status==2, MorphologygroupIndex_corr=="Ductal",  !is.na(ageInt)) |> mutate(outcome=1L)
lcis_c   <- dat |> filter(ethnicityClass==1, status==2, MorphologygroupIndex_corr=="Lobular", !is.na(ageInt)) |> mutate(outcome=1L)
dat_dcis <- bind_rows(ctrl, dcis_c)
dat_lcis <- bind_rows(ctrl, lcis_c)

# ── 4. COHORT SIZE CHECKS ─────────────────────────────────────────────────────
cat("── Cohort sizes ──\n")
chk("BBD cases",             sum(bbd_all$outcome==1, na.rm=TRUE), 1765)
chk("BBD controls",          sum(bbd_all$outcome==0, na.rm=TRUE), 12715)
chk("Non-proliferative BBD", sum(bbd_np$outcome==1,  na.rm=TRUE), 431)
chk("Proliferative BBD",     sum(bbd_pr$outcome==1,  na.rm=TRUE), 274)
chk("DCIS cases",            sum(dat_dcis$outcome==1, na.rm=TRUE), 1663)
chk("DCIS controls",         sum(dat_dcis$outcome==0, na.rm=TRUE), 34424)
chk("LCIS cases",            sum(dat_lcis$outcome==1, na.rm=TRUE), 139)
chk("LCIS controls",         sum(dat_lcis$outcome==0, na.rm=TRUE), 34424)

# Carrier count checks (cross-check against thesis tables)
cat("\n── Carrier counts ──\n")
chk("BRCA2 case carriers / BBD",  sum(bbd_all$BRCA2==1  & bbd_all$outcome==1,  na.rm=TRUE), 11)
chk("ATM   case carriers / DCIS", sum(dat_dcis$ATM==1   & dat_dcis$outcome==1, na.rm=TRUE), 18)
chk("CHEK2 case carriers / DCIS", sum(dat_dcis$CHEK2==1 & dat_dcis$outcome==1, na.rm=TRUE), 40)

# ── 5. TRUNCATING MODELS ──────────────────────────────────────────────────────
cat("\n── Truncating variant models (Firth, Wald CI for speed) ──\n")

r <- fit_or(bbd_all, "BRCA2")
cat(sprintf("  BRCA2/BBD  OR=%.2f p=%.3f  [thesis: 2.87 p=0.015]\n", r$OR, r$p))
chk("BRCA2/BBD OR",  r$OR, 2.87, tol=0.20)
chk("BRCA2/BBD p",   r$p,  0.015, tol=0.008)

r <- fit_or(dat_dcis, "ATM")
cat(sprintf("  ATM/DCIS   OR=%.2f p=%.4f  [thesis: 3.90 p<0.001]\n", r$OR, r$p))
chk("ATM/DCIS OR",   r$OR, 3.90, tol=0.30)
chk("ATM/DCIS p<0.001", as.integer(r$p<0.001), 1L)

r <- fit_or(dat_dcis, "CHEK2")
cat(sprintf("  CHEK2/DCIS OR=%.2f p=%.4f  [thesis: 2.61 p<0.001]\n", r$OR, r$p))
chk("CHEK2/DCIS OR",    r$OR, 2.61, tol=0.20)
chk("CHEK2/DCIS p<0.001", as.integer(r$p<0.001), 1L)

r <- fit_or(dat_dcis, "BRCA2")
cat(sprintf("  BRCA2/DCIS OR=%.2f p=%.3f  [thesis: 2.74 p=0.002]\n", r$OR, r$p))
chk("BRCA2/DCIS OR",  r$OR, 2.74, tol=0.20)
chk("BRCA2/DCIS p",   r$p,  0.002, tol=0.003)

# ── 6. MISSENSE MODELS ────────────────────────────────────────────────────────
cat("\n── Missense variant models ──\n")

r <- fit_or(dat_dcis, "CHEK2_m")
cat(sprintf("  CHEK2_m/DCIS OR=%.2f p=%.4f  [thesis: 2.01 p=0.003]\n", r$OR, r$p))
chk("CHEK2_m/DCIS OR", r$OR, 2.01, tol=0.20)
chk("CHEK2_m/DCIS p",  r$p,  0.003, tol=0.003)

r <- fit_or(dat_dcis, "TP53_m")
cat(sprintf("  TP53_m/DCIS  OR=%.2f p=%.4f  [thesis: 3.65 p=0.006]\n", r$OR, r$p))
chk("TP53_m/DCIS OR",  r$OR, 3.65, tol=0.40)
chk("TP53_m/DCIS p<0.05", as.integer(r$p<0.05), 1L)

# ── 7. FDR CHECK — DCIS truncating ───────────────────────────────────────────
cat("\n── FDR: DCIS truncating (should be FDR-sig: ATM, CHEK2, BRCA2) ──\n")
dcis_p <- map_dfr(GENES, function(g) {
  r <- fit_or(dat_dcis, g)
  tibble(gene=g, OR=r$OR, p=r$p)
}) |> mutate(p_fdr=p.adjust(p,"BH"), fdr_sig=p_fdr<0.05)
sig <- dcis_p |> filter(fdr_sig) |> pull(gene)
cat("  FDR-significant genes:", paste(sort(sig), collapse=", "), "\n")
chk("ATM FDR-sig",   as.integer("ATM"   %in% sig), 1L)
chk("CHEK2 FDR-sig", as.integer("CHEK2" %in% sig), 1L)
chk("BRCA2 FDR-sig", as.integer("BRCA2" %in% sig), 1L)

# ── 8. SUMMARY ────────────────────────────────────────────────────────────────
dir.create("outputs/session_info", showWarnings=FALSE, recursive=TRUE)
si <- sessionInfo()
writeLines(capture.output(si), "outputs/session_info/session_info_verification.txt")

cat(sprintf("\n%s\n", strrep("=",60)))
cat(sprintf("VERIFICATION:  %d PASS  |  %d FAIL\n", PASS, FAIL))
cat(strrep("=",60), "\n")
if (FAIL==0) cat("All checks passed — results fully reproducible.\n") else
  cat("ATTENTION: review failed checks above.\n")
cat(sprintf("R: %s | %s\n", R.version.string, Sys.time()))
