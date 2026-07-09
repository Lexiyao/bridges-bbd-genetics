# Profile-likelihood refit for the two CI-bearing sensitivity figures:
#   (1) ascertainment (5 findings x 3 specs) -> emit profile lo/hi to scratch JSON
#   (2) ATM leave-one-study-out -> update ATM truncating rows in loso_all_fits.rds
# Point estimates are unchanged; only CIs/p move to profile.
# setwd() removed for release — run scripts from the repository root
suppressMessages({library(tidyverse); library(logistf); library(jsonlite)})
GENES <- c("BRCA1","BRCA2","PALB2","CHEK2","ATM","BARD1","RAD51C","RAD51D","TP53")
pheno <- read_delim("concept_807_zhang_bridges_pheno_v17.txt", delim="\t", show_col_types=FALSE,
                    na=c("","NA","888","777","999"))
trunc <- read_csv("concept_807_zhang_bridges_truncating.csv", show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))
miss  <- read_csv("concept_807_zhang_bridges_missense.csv",   show_col_types=FALSE, na=c("","NA")) %>% select(-any_of("...1"))
bc <- function(raw,pat,sfx){cols<-map_chr(GENES,function(g){h<-names(raw)[grepl(paste0("^",g,pat),names(raw),ignore.case=TRUE)];if(length(h)>0)h[1] else NA_character_});k<-!is.na(cols)
  raw%>%select(BRIDGES_ID,all_of(cols[k]))%>%rename(!!!setNames(cols[k],paste0(GENES[k],sfx)))%>%mutate(across(-BRIDGES_ID,~as.integer(replace_na(as.numeric(.x),0)>=1)))}
dat <- pheno %>% filter(ethnicityClass==1,!is.na(ageInt)) %>%
  left_join(bc(trunc,"_truncating$","_t"),by="BRIDGES_ID") %>%
  left_join(bc(miss,"_CADD\\.phred\\.01$","_m"),by="BRIDGES_ID") %>%
  mutate(across(c(ends_with("_t"),ends_with("_m")),~replace_na(.x,0L)))
controls <- dat %>% filter(status==0)
dcis     <- dat %>% filter(status==2, MorphologygroupIndex_corr=="Ductal")
fam <- dat %>% group_by(study) %>% summarise(p=mean(famHist==1,na.rm=TRUE),.groups="drop") %>% filter(p>0.5) %>% pull(study)

pfit <- function(d, col, extra=NULL){
  d<-d%>%mutate(carrier=.data[[col]],study=droplevels(factor(study)))
  if(!is.null(extra)) d<-d%>%filter(!is.na(.data[[extra]]))
  if(sum(d$carrier==1&d$outcome==1)==0) return(c(OR=NA,lo=NA,hi=NA,p=NA))
  form<-if(is.null(extra)) outcome~carrier+ageInt+study else as.formula(paste("outcome~carrier+ageInt+study+",extra))
  f<-tryCatch(logistf(form,data=d,firth=TRUE,pl=TRUE,plconf=2),error=function(e)NULL)
  if(is.null(f)) return(c(OR=NA,lo=NA,hi=NA,p=NA))
  i<-which(names(coef(f))=="carrier")
  c(OR=unname(exp(coef(f)[i])),lo=unname(exp(f$ci.lower[i])),hi=unname(exp(f$ci.upper[i])),p=unname(f$prob[i]))
}
spec_fit <- function(gene,cc,spec){
  ct<-controls; cs<-dcis
  if(spec=="population"){ct<-ct%>%filter(!study%in%fam);cs<-cs%>%filter(!study%in%fam)}
  d<-bind_rows(ct%>%mutate(outcome=0L),cs%>%mutate(outcome=1L))
  pfit(d,paste0(gene,cc),extra=if(spec=="famHist_adjusted")"famHist" else NULL)
}
# ---- (1) ascertainment -> scratch JSON ----
findings <- list(trunc=c("ATM","BRCA2","CHEK2"), miss=c("CHEK2","TP53"))
out <- list(trunc=list(), miss=list())
for (vt in names(findings)) for (g in findings[[vt]]) {
  cc <- if(vt=="miss")"_m" else "_t"
  A<-spec_fit(g,cc,"all_studies"); B<-spec_fit(g,cc,"population"); C<-spec_fit(g,cc,"famHist_adjusted")
  out[[vt]][[g]] <- list(A=A["OR"],A_lo=A["lo"],A_hi=A["hi"],B=B["OR"],B_lo=B["lo"],B_hi=B["hi"],
                         C=C["OR"],C_lo=C["lo"],C_hi=C["hi"])
}
write_json(out, "scratchpad_sens_profile.json", auto_unbox=TRUE, digits=10)
cat("Wrote ascertainment profile JSON.\n")

# ---- (2) ATM LOSO profile -> update rds ----
RDS <- "outputs/tables/loso_all_fits.rds"
if(!file.exists(sub("\\.rds$",".wald.bak.rds",RDS))) file.copy(RDS, sub("\\.rds$",".wald.bak.rds",RDS))
loso <- readRDS(RDS)
case_studies <- sort(unique(dcis$study))
loso_fit <- function(drop=NULL){
  cs<-dcis; ct<-controls
  if(!is.null(drop)){cs<-cs%>%filter(study!=drop); ct<-ct%>%filter(study!=drop)}
  d<-bind_rows(ct%>%mutate(outcome=0L),cs%>%mutate(outcome=1L))%>%
     mutate(carrier=ATM_t, study=droplevels(factor(study)))
  ncc<-sum(d$carrier==1&d$outcome==1); nct<-sum(d$carrier==1&d$outcome==0)
  f<-logistf(outcome~carrier+ageInt+study,data=d,firth=TRUE,pl=TRUE,plconf=2)
  i<-which(names(coef(f))=="carrier")
  tibble(gene="ATM",vtype="truncating",
         dropped=if(is.null(drop))"(none — full data)" else drop,
         n_case=sum(d$outcome==1),carr_case=ncc,carr_ctrl=nct,
         OR=unname(exp(coef(f)[i])),CI_low=unname(exp(f$ci.lower[i])),
         CI_high=unname(exp(f$ci.upper[i])),p_value=unname(f$prob[i]),note="")
}
atm_new <- bind_rows(loso_fit(NULL), map_dfr(case_studies, loso_fit))
loso_other <- loso %>% filter(!(gene=="ATM" & vtype=="truncating"))
# align columns to original
atm_new <- atm_new[, intersect(names(loso), names(atm_new))]
saveRDS(bind_rows(loso_other, atm_new), RDS)
cat("Updated ATM truncating LOSO rows (profile):", nrow(atm_new), "rows\n")
cat("ATM OR range (profile):", sprintf("%.2f-%.2f", min(atm_new$OR[-1]), max(atm_new$OR[-1])), "\n")
