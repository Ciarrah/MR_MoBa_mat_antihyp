`%!in%` = Negate(`%in%`)
library(haven); library(readr)

MBRN=read_sav("")

MBRN=MBRN[,c("PREG_ID_2306","MORS_ALDER","FARS_ALDER", "SVLEN","SVLEN_DG","HYPERTENSJON_KRONISK","HYPERTENSJON_ALENE", "EKLAMPSI","PREEKL","PREEKLTIDL",
             "HELLP", "FSTART","FAAR","FLERFODSEL", "KJONN","VEKT","LENGDE","HODE","APGAR1", "APGAR5")]

# Rename variables

colnames(MBRN)=c("PREG_ID_2306","mthr_age","fthr_age","GA_wks","GA_dys","hyptnsn_pre_prg",""hyptsn_prg","eclmps","preeclmps","early_preeclmps",
                 "HELLP","dlvry_initin","yr_dlvry","sngl_multi_brth","sx_chld","bwt","brth_lngth","head_circum",
                 "apgar_1min","apgar_5min")

colnames(MBRN)[-1]=paste0(colnames(MBRN)[-1],"_MBRN")

lapply(MBRN,summary)
lapply(MBRN,class)

MBRN=lapply(MBRN, as.data.frame)

MBRN[,paste0(c("hyptsn_prg","eclmps","preeclmps","early_preeclmps","HELLP","dlvry_initin"),"_MBRN")]=
  lapply(MBRN[,paste0(c("hyptsn_prg","eclmps","preeclmps","early_preeclmps","HELLP","dlvry_initin"),"_MBRN")], as.factor)

# implausible values e.g. head circumference 3cm w/ bwt 3300g and GA 40 weeks, set implausible values to NA

subset=MBRN[which(MBRN$bwt_MBRN<400&MBRN$dlvry_initin_MBRN%!in%c(2,3)), c("GA_wks_MBRN","bwt_MBRN","brth_lngth_MBRN","head_circum_MBRN","dlvry_initin_MBRN")]
MBRN[which(MBRN$head_circum_MBRN<20&MBRN$GA_wks_MBRN>30),]=NA

# Visually identify outlying points to triangulate

plot(MBRN$bwt_MBRN, MBRN$brth_lngth_MBRN)
plot(MBRN$bwt_MBRN, MBRN$GA_wks_MBRN)
plot(MBRN$GA_wks_MBRN, MBRN$brth_lngth_MBRN)

#1. Set birth length to NA

MBRN[which(MBRN$brth_lngth_MBRN<=10&MBRN$bwt_MBRN>=200),
                 c( "GA_wks_MBRN","bwt_MBRN","brth_lngth_MBRN","head_circum_MBRN","apgar_1min_MBRN")]

MBRN$brth_lngth_MBRN[which(MBRN$brth_lngth_MBRN<=10&MBRN$bwt_MBRN>=200)]=NA

#2. Set birthweight to NA

MBRN[which(MBRN$brth_lngth_MBRN>40&MBRN$bwt_MBRN<=1000),
                 c( "GA_wks_MBRN","bwt_MBRN","brth_lngth_MBRN","head_circum_MBRN","apgar_1min_MBRN")]
MBRN$bwt_MBRN[which(MBRN$brth_lngth_MBRN>40&MBRN$bwt_MBRN<=1000)]=NA
MBRN$bwt_MBRN[which(MBRN$bwt_MBRN==0)]=NA

#3. Set birth length to NA

MBRN$brth_lngth_MBRN[which(MBRN$brth_lngth_MBRN>75)]=NA

#4. Set GA_wks and GA_dys to NA

MBRN[which(MBRN$bwt_MBRN<=1000&MBRN$GA_wks_MBRN>=35),
                 c( "GA_wks_MBRN","bwt_MBRN","brth_lngth_MBRN","head_circum_MBRN","apgar_1min_MBRN")]
MBRN[which(MBRN$bwt_MBRN<=1000&MBRN$GA_wks_MBRN>=35),c("GA_wks_MBRN","GA_dys_MBRN")]=NA

#5. Check plausibility again.

MBRN[which(MBRN$GA_wks_MBRN<=25&MBRN$brth_lngth_MBRN>=40),
                 c( "GA_wks_MBRN","bwt_MBRN","brth_lngth_MBRN","head_circum_MBRN","apgar_1min_MBRN")]

#6. Remove gestational ages greater than 44 weeks - CB implausible.

MBRN[which(MBRN$GA_wks_MBRN>44),c("GA_wks_MBRN","GA_dys_MBRN")]=NA

# Set GA_wks and GA_dys to NA

MBRN[which(MBRN$GA_wks_MBRN<=25&MBRN$brth_lngth_MBRN>=40),c("GA_wks_MBRN","GA_dys_MBRN")]=NA

# Check most implausible outliers removed.

mean(na.omit(MBRN$GA_wks_MBRN))
sd(na.omit(MBRN$GA_wks_MBRN))
range(na.omit(MBRN$GA_wks_MBRN))
summary(MBRN$GA_wks_MBRN)

mean(na.omit(MBRN$bwt_MBRN))
sd(na.omit(MBRN$bwt_MBRN))
range(na.omit(MBRN$bwt_MBRN))
summary(MBRN$bwt_MBRN)

mean(na.omit(MBRN$apgar_1min_MBRN))
sd(na.omit(MBRN$apgar_1min_MBRN))
range(na.omit(MBRN$apgar_1min_MBRN))
summary(MBRN$apgar_1min_MBRN)

mean(na.omit(MBRN$brth_lngth_MBRN))
sd(na.omit(MBRN$brth_lngth_MBRN))
range(na.omit(MBRN$brth_lngth_MBRN))
summary(MBRN$brth_lngth_MBRN)

# Output MBRN summary stats

lapply(MBRN,class)
summary_outcome_2=NULL
summary_outcome_2=as.data.frame(matrix(NA,nrow=(ncol(MBRN)-1),ncol=5))
colnames(summary_outcome_2)=c("variable","mean","sd","proportion","NA")

for (i in colnames(MBRN)[-1])
{
  j=which(colnames(MBRN)[-1]==i)
  if (sapply(MBRN[,i],is.factor)==T)
  {
    summary_outcome_2[j,c("variable","proportion")]=c(i,100*(table(MBRN[,i])/nrow(MBRN)))
    summary_outcome_2[j,"NA"]=length(which(is.na(MBRN[,i])==T))
  }
  else
  {
    summary_outcome_2[j,c("variable","mean","sd")]=c(i,mean(unlist(na.omit(MBRN[,i]))),sd(unlist(na.omit(MBRN[,i]))))
    summary_outcome_2[j,"NA"]=length(which(is.na(MBRN[,i])==T))
  }
}

write.table(summary_outcome_2,"outcome_MBRN_summary.txt")

lapply(MBRN,summary)
lapply(MBRN,class)

MBRN[,c(grep("hyptnsn|HELLP|preeclmps|eclmps|sx_chld", colnames(MBRN)))]=lapply(MBRN[,c(grep("hyptnsn|HELLP|preeclmps|eclmps|sx_chld", colnames(MBRN)))],as.factor)

MBRN$sx_chld_MBRN[which(MBRN$sx_chld_MBRN==0)]=NA
MBRN$sx_chld_MBRN[which(MBRN$sx_chld_MBRN==3)]=NA
MBRN$mthr_age_MBRN[which(MBRN$mthr_age_MBRN>50)]=NA
MBRN$fthr_age_MBRN[which(MBRN$fthr_age_MBRN>60)]=NA

lapply(MBRN,class)

write.csv(MBRN,"MBRN_cleaned.csv")
