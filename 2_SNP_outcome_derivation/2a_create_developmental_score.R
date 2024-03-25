install.packages("misty")
library(haven); library(misty)
`%!in%` = Negate(`%in%`)

# Create the prorated score from item in the six month questionaire. Append to MBRN dataset using child and family ID.

sixmonth=read_sav("")
sixmonth=sixmonth[,c("PREG_ID_2306","BARN_NR","DD949", "DD233", "DD1107", "DD234", "DD1108", "DD950", "DD951")]
sixmonth=as.data.frame(sixmonth)
apply(sixmonth, 2, class)
sixmonth[3:9]=lapply(sixmonth[3:9], factor)

colnames(sixmonth)[3:9]=c("Delayed_development", "Delayed_motor_development", "Delayed_psychomotor_development",
                          "Delayed_motor_development_referred_specialist","Delayed_psychomotor_development_referred_specialist",
                          "Delayed_development_doctor_visit", "Delayed_development_admit_hospital")

apply(sixmonth[3:9], 2, table, useNA = "always")

sixmonth[which(sixmonth$Delayed_motor_development==0),"Delayed_motor_development"]=NA
sixmonth[which(sixmonth$Delayed_motor_development_referred_specialist==0),"Delayed_motor_development_referred_specialist"]=NA
sixmonth$Delayed_motor_development=droplevels(sixmonth$Delayed_motor_development)
sixmonth$Delayed_motor_development_referred_specialist=droplevels(sixmonth$Delayed_motor_development_referred_specialist)

levels(sixmonth$Delayed_development)=levels(sixmonth$Delayed_motor_development)=levels(sixmonth$Delayed_psychomotor_development)=
  levels(sixmonth$Delayed_development_doctor_visit)=levels(sixmonth$Delayed_development_admit_hospital)=
  c("No", "Yes")
levels(sixmonth$Delayed_motor_development_referred_specialist)=levels(sixmonth$Delayed_psychomotor_development_referred_specialist)=
  c("No", "Yes, referred from health centre", "Yes, referred by someone else")

apply(sixmonth[,3:9], 2, table, useNA = "always")

sixmonth$any_development_delay=ifelse(sixmonth$Delayed_development =="Yes" | sixmonth$Delayed_motor_development =="Yes" |
                                        sixmonth$Delayed_psychomotor_development =="Yes", "Yes", NA)

sixmonth[which(is.na(sixmonth$any_development_delay)&(sixmonth$Delayed_development =="No" | sixmonth$Delayed_motor_development =="No" |
                                                        sixmonth$Delayed_psychomotor_development =="No")),"any_development_delay"]="No"
sixmonth=sixmonth[,c("PREG_ID_2306", "BARN_NR", "any_development_delay")]

colnames(MBRN_families)[grep("BARN",colnames(MBRN_families))]=gsub("_MBRN","",colnames(MBRN_families[grep("BARN",colnames(MBRN_families))]))
outcomes=merge(MBRN_families, sixmonth, by=c("PREG_ID_2306","BARN_NR"), all=T)
rm(MBRN_families)

sixmonth=read_sav("")
sixmonth=sixmonth[,c("PREG_ID_2306","BARN_NR","DD348", "DD349", "DD350", "DD351", "DD352", "DD353", "DD354", "DD355", "DD356",
                     "DD357", "DD358")]

# replace 1 with 10, 2 with 5, 3 with 0 and 4 with NA

rplc_w_tn=colnames(sixmonth)[grep("DD", colnames(sixmonth))]
for (i in rplc_w_tn)
{
  sixmonth[which(sixmonth[,i]==1),i]=10
  sixmonth[which(sixmonth[,i]==2),i]=5
  sixmonth[which(sixmonth[,i]==3),i]=0
  sixmonth[which(sixmonth[,i]%in%c(0,4)),i]=NA
}

# Prorated score calculation - calculate the number of NA values in each row, max of 3 allowed (AH). Otherwise set NA.

non_na_counts=rowSums(!is.na(sixmonth[,3:13]))
sixmonth$prorated_scores=NA
for (i in 1:nrow(sixmonth)) {
  if (non_na_counts[i] >= 3) {
    sixmonth[i,"prorated_scores"]=item.scores(sixmonth[i,3:13])
  } else {
    sixmonth[i,"prorated_scores"]=NA
  }
}

outcomes=merge(outcomes, sixmonth[,c("PREG_ID_2306", "BARN_NR", "prorated_scores")], all.x=T)

write.csv(outcomes, "MBRN_sixmonth_cleaned.csv")
