# =============================================================================
# Unify CI method on profile-likelihood (Path 1)
# Re-runs the truncating primary table (4.2) and the familial/population
# sensitivity table (4.5) with logistf pl = TRUE, to match the missense
# primary table (4.3, already pl = TRUE). Prints Wald vs profile side-by-side
# so the point estimates can be confirmed identical (only CI/p move).
#
# plconf = carrier index only -> profiles just the carrier coefficient.
# This is IDENTICAL to the original missense run's carrier CI (logistf default
# profiles all coefs; the carrier CI does not depend on plconf) but far faster.
# =============================================================================
# setwd() removed for release — run scripts from the repository root
suppressMessages({library(tidyverse); library(logistf)})
options(width = 200)

PHENO_FILE <- "concept_807_zhang_bridges_pheno_v17.txt"
TRUNC_FILE <- "concept_807_zhang_bridges_truncating.csv"
MISS_FILE  <- "concept_807_zhang_bridges_missense.csv"
GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")

pheno <- read_delim(PHENO_FILE, delim="\t", show_col_types=FALSE,
                    na=c("","NA","888","777","999"))
trunc <- read_csv(TRUNC_FILE, show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))
miss  <- read_csv(MISS_FILE,  show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))

build_carriers <- function(raw, pattern, suffix) {
  cols <- map_chr(GENES, function(g) {
    hit <- names(raw)[grepl(paste0("^", g, pattern), names(raw), ignore.case=TRUE)]
    if (length(hit)>0) hit[1] else NA_character_ })
  keep <- !is.na(cols)
  raw %>% select(BRIDGES_ID, all_of(cols[keep])) %>%
    rename(!!!setNames(cols[keep], paste0(GENES[keep], suffix))) %>%
    mutate(across(-BRIDGES_ID, ~ as.integer(replace_na(as.numeric(.x),0) >= 1)))
}
carriers_t <- build_carriers(trunc, "_truncating$",      "_t")
carriers_m <- build_carriers(miss,  "_CADD\\.phred\\.01$","_m")

dat <- pheno %>% filter(ethnicityClass==1, !is.na(ageInt)) %>%
  left_join(carriers_t, by="BRIDGES_ID") %>%
  left_join(carriers_m, by="BRIDGES_ID") %>%
  mutate(across(c(ends_with("_t"), ends_with("_m")), ~ replace_na(.x, 0L)))

# Firth fit returning BOTH Wald and profile CI for the carrier term -----------
fit_both <- function(d, carrier_col, extra = NULL) {
  d <- d %>% mutate(carrier = .data[[carrier_col]], study = factor(study))
  ncc <- sum(d$carrier==1 & d$outcome==1); nct <- sum(d$carrier==1 & d$outcome==0)
  if (ncc + nct == 0)
    return(tibble(carr_case=ncc, carr_ctrl=nct, OR=NA, w_lo=NA,w_hi=NA,w_p=NA,
                  p_lo=NA,p_hi=NA,p_p=NA))
  form <- if (is.null(extra)) outcome ~ carrier + ageInt + study
          else as.formula(paste("outcome ~ carrier + ageInt + study +", extra))
  if (!is.null(extra)) d <- d %>% filter(!is.na(.data[[extra]]))
  fw <- logistf(form, data=d, firth=TRUE, pl=FALSE)
  fp <- logistf(form, data=d, firth=TRUE, pl=TRUE, plconf=2)  # carrier = idx 2
  i <- which(names(coef(fw))=="carrier")
  tibble(carr_case=ncc, carr_ctrl=nct,
         OR   = exp(coef(fw)[i]),
         w_lo = exp(fw$ci.lower[i]), w_hi = exp(fw$ci.upper[i]), w_p = fw$prob[i],
         p_lo = exp(fp$ci.lower[i]), p_hi = exp(fp$ci.upper[i]), p_p = fp$prob[i])
}

panel <- function(d, cc, lab) map_dfr(GENES, function(g)
  fit_both(d, paste0(g, cc)) %>% mutate(gene=g, analysis=lab)) %>%
  mutate(w_fdr = p.adjust(w_p, "BH"), p_fdr = p.adjust(p_p, "BH"))

# Cohorts ---------------------------------------------------------------------
controls <- dat %>% filter(status==0) %>% mutate(outcome=0L)
dcis     <- dat %>% filter(status==2, MorphologygroupIndex_corr=="Ductal")  %>% mutate(outcome=1L)
lcis     <- dat %>% filter(status==2, MorphologygroupIndex_corr=="Lobular") %>% mutate(outcome=1L)
bbd_all  <- dat %>% filter(status==0, !is.na(BBD_history)) %>% mutate(outcome=as.integer(BBD_history==1))

dat_dcis <- bind_rows(controls, dcis)
dat_lcis <- bind_rows(controls, lcis)

show <- function(df, title) {
  cat("\n", strrep("=",110), "\n", title, "\n", strrep("=",110), "\n", sep="")
  df %>% transmute(gene,
      `carr(case/ctrl)` = sprintf("%d/%d", carr_case, carr_ctrl),
      OR = ifelse(is.na(OR),NA,sprintf("%.2f",OR)),
      `Wald CI`    = ifelse(is.na(OR),"—",sprintf("%.2f-%.2f", w_lo, w_hi)),
      `Profile CI` = ifelse(is.na(OR),"—",sprintf("%.2f-%.2f", p_lo, p_hi)),
      `Wald q`    = sprintf("%.3f", w_fdr),
      `Profile q` = sprintf("%.3f", p_fdr)) %>% as.data.frame() %>% print(row.names=FALSE)
}

cat("\n############ TABLE 4.2  TRUNCATING — Wald (current) vs Profile (unified) ############\n")
show(panel(bbd_all,  "_t", "BBD"),  "Table 4.2  BBD  (truncating)")
show(panel(dat_dcis, "_t", "DCIS"), "Table 4.2  DCIS (truncating)")
show(panel(dat_lcis, "_t", "LCIS"), "Table 4.2  LCIS (truncating)")

cat("\n############ TABLE 4.3  MISSENSE DCIS — confirm already-profile is reproduced ############\n")
show(panel(dat_dcis, "_m", "DCIS miss"), "Table 4.3  DCIS (missense)")

# Table 4.5 sensitivity: 5 findings x 3 specs --------------------------------
fam_studies <- dat %>% group_by(study) %>%
  summarise(pct=mean(famHist==1,na.rm=TRUE),.groups="drop") %>%
  filter(pct>0.50) %>% pull(study)
cat("\nFamilial-ascertainment studies:", paste(fam_studies, collapse=", "), "\n")

sens <- function(gene, cc, spec) {
  ct <- controls; cs <- dcis
  if (spec=="population") { ct<-ct%>%filter(!study%in%fam_studies); cs<-cs%>%filter(!study%in%fam_studies) }
  d <- bind_rows(ct%>%mutate(outcome=0L), cs%>%mutate(outcome=1L))
  fit_both(d, paste0(gene,cc), extra = if(spec=="famHist_adjusted") "famHist" else NULL) %>%
    mutate(finding=paste(gene,cc), spec=spec)
}
findings <- list(c("ATM","_t"),c("BRCA2","_t"),c("CHEK2","_t"),c("CHEK2","_m"),c("TP53","_m"))
cat("\n############ TABLE 4.5  SENSITIVITY — Wald (current) vs Profile (unified) ############\n")
for (f in findings) for (s in c("all_studies","population","famHist_adjusted")) {
  r <- sens(f[1], f[2], s)
  cat(sprintf("%-12s %-17s OR=%5.2f | Wald %5.2f-%6.2f | Profile %5.2f-%6.2f\n",
              r$finding, s, r$OR, r$w_lo, r$w_hi, r$p_lo, r$p_hi))
}
cat("\nDONE.\n")
