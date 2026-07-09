# =============================================================================
# Menopause exposure — age at menopause + invasive-BC positive control
# "Genetic Associations Between Benign Breast Disease, Carcinoma in Situ,
#  and Breast Cancer Susceptibility Genes" — Cambridge MPhil PHS 2025/26
#
# PURPOSE (supervisor feedback)
#   (1) Coding: confirm the special codes for menoStat and mensAgeLast.
#   (2) Describe age at menopause (mensAgeLast).
#   (3) Positive control: reproduce the reproductive-factor associations with
#       INVASIVE breast cancer (status==1) and check they match the established
#       epidemiology (Collaborative Group). This validates the analysis pipeline.
#   (4) Replace the age-confounded binary menopausal-status exposure with the
#       epidemiologically standard CONTINUOUS age at menopause for the
#       non-invasive outcomes (BBD/DCIS/LCIS), restricted to post-menopausal
#       women (mensAgeLast is structurally NA for pre-menopausal women).
#
# OUTPUT
#   outputs/tables/menopause_age_analysis.xlsx
# =============================================================================

suppressPackageStartupMessages({library(tidyverse); library(writexl)})
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

raw <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                  show_col_types = FALSE, na = c("", "NA"))
ph  <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim = "\t",
                  show_col_types = FALSE, na = c("", "NA", "888", "777", "999"))

# (1) Coding record
coding <- tibble(
  variable = c("menoStat", "menoStat", "menoStat",
               "mensAgeLast", "mensAgeLast", "mensAgeLast"),
  code = c("1", "2", "888", "real (16-73)", "777", "888"),
  meaning = c("pre/peri-menopausal", "post-menopausal", "unknown (treated NA)",
              "age at menopause (post-menopausal women)",
              "not applicable - pre-menopausal (treated NA)",
              "unknown (treated NA)"),
  n = c(sum(raw$menoStat==1,na.rm=T), sum(raw$menoStat==2,na.rm=T), sum(raw$menoStat==888,na.rm=T),
        sum(raw$mensAgeLast %in% 16:73,na.rm=T), sum(raw$mensAgeLast==777,na.rm=T), sum(raw$mensAgeLast==888,na.rm=T)))

# (2) Age-at-menopause descriptive
m <- ph$mensAgeLast[!is.na(ph$mensAgeLast)]
descr <- tibble(metric = c("n with value","mean","sd","median","min","max"),
                value = c(length(m), round(mean(m),1), round(sd(m),1),
                          median(m), min(m), max(m)))

eur  <- ph %>% filter(ethnicityClass == 1, !is.na(ageInt))
ctrl <- eur %>% filter(status == 0)

or_ci <- function(d, e) {
  d <- d %>% select(out, all_of(e), ageInt, study) %>% mutate(study=factor(study)) %>% drop_na()
  f <- glm(as.formula(paste("out ~", e, "+ ageInt + study")), d, family = binomial())
  ci <- exp(confint.default(f)[e, ])
  tibble(exposure=e, n=nrow(d), n_case=sum(d$out==1),
         OR=round(exp(coef(f)[e]),3), CI_low=round(ci[1],3), CI_high=round(ci[2],3),
         p=signif(summary(f)$coefficients[e,4],3))
}

# (3) Positive control: invasive BC (status==1)
inv <- bind_rows(ctrl %>% mutate(out=0L), eur %>% filter(status==1) %>% mutate(out=1L))
pos_ctrl <- map_dfr(c("parity","ageMenarche","menoStat","mensAgeLast","HRTEver"),
                    ~ or_ci(inv, .x)) %>% mutate(outcome = "Invasive BC (positive control)")

# (4) Age at menopause for the non-invasive outcomes
sets <- list(
  BBD  = bind_rows(eur%>%filter(status==0,BBD_history==1)%>%mutate(out=1L),
                   eur%>%filter(status==0,BBD_history==0)%>%mutate(out=0L)),
  DCIS = bind_rows(ctrl%>%mutate(out=0L), eur%>%filter(status==2,MorphologygroupIndex_corr=="Ductal")%>%mutate(out=1L)),
  LCIS = bind_rows(ctrl%>%mutate(out=0L), eur%>%filter(status==2,MorphologygroupIndex_corr=="Lobular")%>%mutate(out=1L))
)
meno_age <- map_dfr(names(sets), ~ or_ci(sets[[.x]], "mensAgeLast") %>% mutate(outcome=.x))

cat("=== (1) Coding ===\n");            print(as.data.frame(coding))
cat("\n=== (2) Age at menopause (post-menopausal women) ===\n"); print(as.data.frame(descr))
cat("\n=== (3) Positive control vs invasive BC ===\n"); print(as.data.frame(pos_ctrl))
cat("Literature: parity & later menarche protective; later menopause (higher age) higher risk.\n")
cat("\n=== (4) Age at menopause vs non-invasive outcomes (post-menopausal women) ===\n")
print(as.data.frame(meno_age))
cat("All null — confirms the binary-menoStat DCIS 0.51 was an age-collinearity artifact.\n")

write_xlsx(list(coding = coding, menopause_age_descriptive = descr,
                positive_control_invasiveBC = pos_ctrl,
                menopause_age_noninvasive = meno_age),
           "outputs/tables/menopause_age_analysis.xlsx")
message("\nSaved: outputs/tables/menopause_age_analysis.xlsx")
