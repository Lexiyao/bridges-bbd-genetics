# Regenerate DCIS_LCIS_truncating_results.xlsx with profile-likelihood CIs (pl=TRUE).
# Backs up the original (Wald) workbook first. Schema identical to the original so
# build_exhibits.R reads it unchanged.
# setwd() removed for release — run scripts from the repository root
suppressMessages({library(tidyverse); library(logistf); library(writexl)})
OUT <- "outputs/tables/DCIS_LCIS_truncating_results.xlsx"
if (file.exists(OUT) && !file.exists(sub("\\.xlsx$",".wald.bak.xlsx",OUT)))
  file.copy(OUT, sub("\\.xlsx$",".wald.bak.xlsx",OUT))

GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")
pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim="\t", show_col_types=FALSE,
                    na=c("","NA","888","777","999"))
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv", show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))
cols <- map_chr(GENES,function(g){h<-names(trunc)[grepl(paste0("^",g,"_truncating$"),names(trunc),ignore.case=TRUE)];if(length(h)>0)h[1] else NA_character_})
tsel <- trunc %>% select(BRIDGES_ID, all_of(cols)) %>% rename(!!!setNames(cols,GENES)) %>%
  mutate(across(-BRIDGES_ID, ~as.integer(replace_na(as.numeric(.x),0)>=1)))
dat <- pheno %>% left_join(tsel,by="BRIDGES_ID") %>% mutate(across(all_of(GENES),~replace_na(.x,0L)))

controls <- dat %>% filter(status==0, ethnicityClass==1, !is.na(ageInt)) %>% mutate(outcome=0L)
dcis <- dat %>% filter(status==2, ethnicityClass==1, MorphologygroupIndex_corr=="Ductal", !is.na(ageInt)) %>% mutate(outcome=1L)
lcis <- dat %>% filter(status==2, ethnicityClass==1, MorphologygroupIndex_corr=="Lobular", !is.na(ageInt)) %>% mutate(outcome=1L)

run_firth <- function(data, gene){
  ncc<-sum(data[[gene]]==1 & data$outcome==1); nct<-sum(data[[gene]]==1 & data$outcome==0)
  nca<-sum(data$outcome==1); ncn<-sum(data$outcome==0)
  if(ncc+nct==0) return(tibble(gene,n_case=nca,n_ctrl=ncn,n_carrier_case=0L,n_carrier_ctrl=0L,
                               OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p_value=NA_real_,method="no carriers"))
  md <- data %>% select(outcome, carrier=all_of(gene), ageInt, study) %>% mutate(study=factor(study)) %>% drop_na()
  fit <- tryCatch(logistf(outcome~carrier+ageInt+study, data=md, firth=TRUE, pl=TRUE, plconf=2), error=function(e)NULL)
  if(is.null(fit)) return(tibble(gene,n_case=nca,n_ctrl=ncn,n_carrier_case=ncc,n_carrier_ctrl=nct,
                                 OR=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p_value=NA_real_,method="failed"))
  i<-which(names(coef(fit))=="carrier")
  tibble(gene,n_case=nca,n_ctrl=ncn,n_carrier_case=ncc,n_carrier_ctrl=nct,
         OR=unname(exp(coef(fit)[i])), CI_low=unname(exp(fit$ci.lower[i])), CI_high=unname(exp(fit$ci.upper[i])),
         p_value=unname(fit$prob[i]), method="Firth")
}
panel <- function(data,label) map_dfr(GENES,~run_firth(data,.x)) %>%
  mutate(analysis=label, p_fdr=p.adjust(p_value,"BH"), fdr_sig=p_fdr<0.05) %>%
  select(analysis,gene,n_case,n_ctrl,n_carrier_case,n_carrier_ctrl,OR,CI_low,CI_high,p_value,p_fdr,fdr_sig,method)

res_dcis <- panel(bind_rows(controls,dcis), "DCIS vs Controls (truncating)")
res_lcis <- panel(bind_rows(controls,lcis), "LCIS vs Controls (truncating)")

fmt <- function(df) df %>% transmute(Gene=gene,`N (cases)`=n_case,`N (controls)`=n_ctrl,
  `Carriers (cases)`=n_carrier_case,`Carriers (controls)`=n_carrier_ctrl,
  `OR (95% CI)`=ifelse(is.na(OR),"—",sprintf("%.2f (%.2f–%.2f)",OR,CI_low,CI_high)),
  `p-value`=ifelse(is.na(p_value),"—",ifelse(p_value<0.001,"<0.001",sprintf("%.3f",p_value))),
  `p-FDR (BH)`=ifelse(is.na(p_fdr),"—",ifelse(p_fdr<0.001,"<0.001",sprintf("%.3f",p_fdr))),
  `FDR sig.`=ifelse(is.na(fdr_sig),"—",ifelse(fdr_sig,"Yes*","No")))

write_xlsx(list(DCIS_truncating=fmt(res_dcis), LCIS_truncating=fmt(res_lcis),
                Raw_results=bind_rows(res_dcis,res_lcis)), path=OUT)
cat("Rewrote (profile):",OUT,"\n")
print(as.data.frame(bind_rows(res_dcis,res_lcis) %>% select(analysis,gene,OR,CI_low,CI_high,fdr_sig)), row.names=FALSE)
