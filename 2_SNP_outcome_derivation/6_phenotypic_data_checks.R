`%!in%` = Negate(`%in%`)
library(data.table) 

pheno_geno=fread("pheno_geno_pca_merged.txt")
pheno_geno=as.data.frame(pheno_geno)

# Missing data fix - set variables that are binary to be NA if not 1.

tmp=grep("hyptnsn_pre_prg|plnd_c_sctn|prvs_c_sctn|prgncy_cmplc|imp_ftl_grwth|prfrnc|c_sctn|brch|othr|brth_cmplctn|
         shyptsn|eclmps|HELLP|hyptsn",colnames(pheno_geno),value=T)

for (i in tmp)
{
  pheno_geno[which(is.na(pheno_geno[,i])),i]=0
}

factors=c("hyptnsn_pre_prg_MBRN","eplps_MBRN", "hyptsn_prg_MBRN","eclmps_MBRN","preeclmps_MBRN", "early_preeclmps_MBRN","HELLP_MBRN",
          "sx_chld_MBRN","dlvry_initin_MBRN")

numeric=c("head_circum_MBRN","apgar_1min_MBRN","apgar_5min_MBRN", "bwt_MBRN", "brth_lngth_MBRN")

# Binary variables: number of NA, proportion and counts for each level.

for(j in factors) set(pheno_geno, j = j, value = factor(pheno_geno[[j]], exclude = NULL))

sink(file("binary_summary_stats.txt"), type="output")
for (i in 1:length(factors))
{
  print(factors[i])
  print(summary(pheno_geno[,factors[i]]))
  print(round((summary(pheno_geno[,factors[i]])/sum(summary(pheno_geno[,factors[i]])))*100,2))
}
sink()

# Continuous variable: mean, sd and count of NA.

for (j in numeric) set(pheno_geno, j = j, value = as.numeric(pheno_geno[[j]], exclude = NULL))

sink(file("continuous_summary_stats.txt"), type="output")
for (i in 1:length(numeric))
{
  print(numeric[i])
  print(c("mean:",round(mean(na.omit(pheno_geno[,numeric[i]])),2)))
  print(c("sd:",round(sd(na.omit(pheno_geno[,numeric[i]])),2)))
  print(c("is na:",table(is.na(pheno_geno[,numeric[i]]))))
}
sink()

fwrite(pheno_geno,"pheno_geno_hyp_post_summary_stats.txt")
