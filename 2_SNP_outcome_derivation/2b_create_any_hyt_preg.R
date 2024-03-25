library(data.table)

pheno_geno=fread("pheno_geno_hyp.txt")

outcome_list=grep("hyptsn_prg_MBRN|eclmps_MBRN|preeclmps_MBRN|early_preeclmps_MBRN|HELLP_MBRN",colnames(pheno_geno), value=T)
pheno_geno$any_hyp_prg=NA
pheno_geno$any_hyp_prg=ifelse(rowSums(as.data.frame(pheno_geno)[,paste0(outcome_list)]!=0)>0,1,0)

table(pheno_geno$any_hyp_prg)
pheno_geno$any_hyp_prg=as.factor(pheno_geno$any_hyp_prg)

outcomes=read.csv("MBRN_sixmonth_cleaned.csv")
outcomes=merge(outcomes, pheno_geno[,c("PREG_ID_2306", "BARN_NR", "any_hyp_prg")], all.x=T)
