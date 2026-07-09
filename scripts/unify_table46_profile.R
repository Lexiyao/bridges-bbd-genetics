# Recompute Table 4.6 (multiple-testing sensitivity) under profile-likelihood.
# Profile (penalised-LRT) p-values for all 54 tests (9 genes x 2 variant x 3 outcome),
# then per-family BH (within outcome x variant, 9 each), global BH (54), Bonferroni (x54).
# setwd() removed for release — run scripts from the repository root
suppressMessages({library(tidyverse); library(logistf)})
options(width = 200)
GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")

pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim="\t",
                    show_col_types=FALSE, na=c("","NA","888","777","999"))
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv", show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))
miss  <- read_csv("concept_807_zhang_bridges_missense.csv",   show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))
bc <- function(raw, pat, sfx){cols<-map_chr(GENES,function(g){h<-names(raw)[grepl(paste0("^",g,pat),names(raw),ignore.case=TRUE)];if(length(h)>0)h[1] else NA_character_});k<-!is.na(cols)
  raw%>%select(BRIDGES_ID,all_of(cols[k]))%>%rename(!!!setNames(cols[k],paste0(GENES[k],sfx)))%>%mutate(across(-BRIDGES_ID,~as.integer(replace_na(as.numeric(.x),0)>=1)))}
dat <- pheno %>% filter(ethnicityClass==1,!is.na(ageInt)) %>%
  left_join(bc(trunc,"_truncating$","_t"),by="BRIDGES_ID") %>%
  left_join(bc(miss,"_CADD\\.phred\\.01$","_m"),by="BRIDGES_ID") %>%
  mutate(across(c(ends_with("_t"),ends_with("_m")),~replace_na(.x,0L)))

controls<-dat%>%filter(status==0)%>%mutate(outcome=0L)
dcis<-dat%>%filter(status==2,MorphologygroupIndex_corr=="Ductal")%>%mutate(outcome=1L)
lcis<-dat%>%filter(status==2,MorphologygroupIndex_corr=="Lobular")%>%mutate(outcome=1L)
bbd<-dat%>%filter(status==0,!is.na(BBD_history))%>%mutate(outcome=as.integer(BBD_history==1))
cohorts<-list(BBD=bind_rows(controls%>%filter(!is.na(BBD_history))%>%mutate(outcome=0L)%>%filter(BBD_history!=1|is.na(BBD_history)),NULL))
# simpler: BBD cohort already has outcome 0/1 within bbd
pl_p <- function(d, col){
  d<-d%>%mutate(carrier=.data[[col]],study=factor(study))
  ncc<-sum(d$carrier==1&d$outcome==1); nct<-sum(d$carrier==1&d$outcome==0)
  if(ncc+nct==0) return(c(OR=NA,p=NA))
  f<-tryCatch(logistf(outcome~carrier+ageInt+study,data=d,firth=TRUE,pl=TRUE,plconf=2),error=function(e)NULL)
  if(is.null(f)) return(c(OR=NA,p=NA))
  i<-which(names(coef(f))=="carrier"); c(OR=unname(exp(coef(f)[i])),p=unname(f$prob[i]))
}
panels <- list(
  list(o="BBD", v="truncating", d=bbd,                      cc="_t"),
  list(o="DCIS",v="truncating", d=bind_rows(controls,dcis), cc="_t"),
  list(o="LCIS",v="truncating", d=bind_rows(controls,lcis), cc="_t"),
  list(o="BBD", v="missense",   d=bbd,                      cc="_m"),
  list(o="DCIS",v="missense",   d=bind_rows(controls,dcis), cc="_m"),
  list(o="LCIS",v="missense",   d=bind_rows(controls,lcis), cc="_m"))
res <- map_dfr(panels, function(P) map_dfr(GENES, function(g){
  r<-pl_p(P$d, paste0(g,P$cc)); tibble(outcome=P$o,variant=P$v,gene=g,OR=unname(r["OR"]),p=unname(r["p"]))}))

res <- res %>% group_by(outcome,variant) %>% mutate(q_perfam=p.adjust(p,"BH")) %>% ungroup() %>%
  mutate(q_global=p.adjust(p,"BH"), p_bonf=pmin(p*sum(!is.na(p)),1))
cat("Total tests with a p-value:", sum(!is.na(res$p)), "\n\n")

findings <- tribble(~outcome,~variant,~gene,
  "DCIS","truncating","ATM","DCIS","truncating","BRCA2","DCIS","truncating","CHEK2",
  "DCIS","missense","CHEK2","DCIS","missense","TP53")
tab46 <- findings %>% left_join(res,by=c("outcome","variant","gene")) %>%
  mutate(survives=case_when(p_bonf<0.05~"Bonferroni", q_global<0.05~"global BH",
                            q_perfam<0.05~"per-family BH", TRUE~"none")) %>%
  transmute(Finding=paste0(gene," (",ifelse(variant=="truncating","PTV","missense"),")"),
            OR=sprintf("%.2f",OR), p=sprintf("%.4f",p),
            `per-family BH q`=sprintf("%.3f",q_perfam),
            `global BH q (54)`=sprintf("%.3f",q_global),
            `Bonferroni p`=sprintf("%.3f",p_bonf), `Survives up to`=survives)
cat("############ TABLE 4.6  UNIFIED (profile-likelihood) ############\n")
print(as.data.frame(tab46), row.names=FALSE)
cat("\n--- old (mixed Wald/profile) for comparison ---\n")
cat("ATM PTV   : perfam 0.000 | global 0.000 | bonf 0.000 | Bonferroni\n")
cat("BRCA2 PTV : perfam 0.004 | global 0.025 | bonf 0.074 | global BH\n")
cat("CHEK2 PTV : perfam 0.001 | global 0.007 | bonf 0.015 | Bonferroni\n")
cat("CHEK2 miss: perfam 0.023 | global 0.034 | bonf 0.137 | global BH\n")
cat("TP53 miss : perfam 0.027 | global 0.065 | bonf 0.324 | per-family BH\n")
