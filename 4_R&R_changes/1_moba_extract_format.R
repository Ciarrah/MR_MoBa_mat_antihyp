library(tidyverse);library(stringr);library(readr);library(haven);library(data.table);library(misty)
`%!in%` = Negate(`%in%`)

MBRN=read_sav("")
MBRN=MBRN[,c("PREG_ID_2306","BARN_NR","MORS_ALDER","FARS_ALDER", "ZSCORE_BW_GA")]

colnames(MBRN)=c("PREG_ID_2306","BARN_NR","mthr_age","fthr_age", "BW_GA",)

# Export to cluster to format

write.csv(MBRN,"r_r_phenotypic_data.csv")

rm(MBRN)

# Obtain merged file from cluster.

genetic_data=read_csv("trio_subset_hyp.csv")

# Load in principal components to merge in

pca=read_table2("MoBaPsychGen_v1-ec-eur-batch-basic-qc-cov.txt",
                col_names = T)
pca$ID_2306[which(pca$Role=="Father")]=gsub("F","",pca$ID_2306[which(pca$Role=="Father")])
pca$ID_2306[which(pca$Role=="Mother")]=gsub("M","",pca$ID_2306[which(pca$Role=="Mother")])

# Split id_2306 into BARN.NR (birth registry file) and preg ID and merge both

child_details=as.data.frame(matrix(NA,nrow=length((pca$ID_2306[which(pca$Role=="Child")])),ncol=3))
child_details$V1=unlist(map(str_split(pca$ID_2306[which(pca$Role=="Child")],"_"),1))
child_details$V2=unlist(map(str_split(pca$ID_2306[which(pca$Role=="Child")],"_"),2))
child_details$V3=pca$IID[which(pca$Role=="Child")]
colnames(child_details)=c("ID_2306","BARN_NR","IID")

# Duplicated genotyping occurred for random people, of which a random one was selected and included so there 
# are only unique member in the sample. Multiple genotyped people have a suffix e.g. _D1 etc.

FID=pca[which(pca$Role=="Father"),c("ID_2306","IID")]
MID=pca[which(pca$Role=="Mother"),c("ID_2306","IID")]
Linkage_merge_file=merge(merge(FID,MID,by="ID_2306",all=T),child_details,by="ID_2306",all=T)
rm(FID,MID,child_details)
colnames(Linkage_merge_file)=c("PREG_ID_2306","PAT","MAT","BARN_ID","IID")
names(pca)[names(pca) == 'ID_2306']='PREG_ID_2306'

phenotypic=read.csv("r_r_phenotypic_data.csv")
pheno_w_link_merge=merge(phenotypic,Linkage_merge_file,by="PREG_ID_2306")
rm(phenotypic, Linkage_merge_file)
pheno_w_link_merge[,c("MAT","PAT")]=NULL

genetic_data[,c("MAT","PAT")]=NULL
geno_pheno_w_link_merge=merge(pheno_w_link_merge, genetic_data, by="IID", all=T)
rm(pheno_w_link_merge, genetic_data)
check=apply(geno_pheno_w_link_merge[,grep("_offspring|_mother|_father",colnames(geno_pheno_w_link_merge))], 1, function(x)all(is.na(x)))
check_2=apply(geno_pheno_w_link_merge[,grep("_offspring|_mother|_father",colnames(geno_pheno_w_link_merge))], 1, function(x)any(is.na(x)))
geno_pheno_w_link_merge=geno_pheno_w_link_merge[which(check==F),]
pca_geno_pheno_w_link_merge=merge(geno_pheno_w_link_merge,pca,by="IID")
pca_geno_pheno_w_link_merge=distinct(pca_geno_pheno_w_link_merge)
write.table(pca_geno_pheno_w_link_merge,"pheno_geno_pca_merged_r_r_hyp.txt")

# Summary stats

pheno_geno=fread("pheno_geno_pca_merged_r_r_hyp.txt")
pheno_geno=as.data.frame(pheno_geno)

numeric="BW_GA"

# continuous variable need mean, sd and count of NA

for (j in numeric) set(pheno_geno, j = j, value = as.numeric(pheno_geno[[j]], exclude = NULL))

sink(file("continuous_summary_stats_r_r.txt"), type="output")
for (i in 1:length(numeric))
{
  print(numeric[i])
  print(c("mean:",round(mean(na.omit(pheno_geno[,numeric[i]])),2)))
  print(c("sd:",round(sd(na.omit(pheno_geno[,numeric[i]])),2)))
  print(c("is na:",table(is.na(pheno_geno[,numeric[i]]))))
}
sink()

fwrite(pheno_geno,"pheno_geno_r_r_summary_stats.txt")

## Regression

outcome_list="BW_GA"

snps=c("rs10764331","rs12258967","rs1262894","rs13143677","rs1557765","rs1801253","rs3821843")

lm_MAT=lm_PAT=lm_OFF=as.data.frame(matrix(NA,nrow=7,ncol=4))
lm_ALL=as.data.frame(matrix(NA,nrow=7,ncol=10))
MAT=PAT=OFF=ALL=NULL
colnames(lm_MAT)=colnames(lm_PAT)=colnames(lm_OFF)=paste0(c("Beta","SE","Pval","DF"),"_")
              
colnames(lm_ALL)=paste0(c(paste0(c(rep("Beta",3),rep("SE",3),rep("Pval",3)),rep(c("_Off","_Mat","_Pat"),3)),"DF"),"_",
                        c(rep("BW_GA",each=10)))

rownames(lm_MAT)=rownames(lm_PAT)=rownames(lm_OFF)=rownames(lm_ALL)=snps

index_1=1
index_2=1

for (k in snps)
{
  MAT=grep(k,grep("mothers",colnames(pheno_geno),value=T),value=T)
  PAT=grep(k,grep("fathers",colnames(pheno_geno),value=T),value=T)
  OFF=grep(k,grep("offspring",colnames(pheno_geno),value=T),value=T)
  for (j in 1:length(outcome_list))
  {
    if ((is.factor(pheno_geno[[outcome_list[j]]]))==T&&length(levels(pheno_geno[[outcome_list[j]]]))==2)
    {
      lm_MAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(MAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                              sep="~")),data=pheno_geno,family="binomial"))
      lm_PAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                              sep="~")),data=pheno_geno,family="binomial"))
      lm_OFF_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                              sep="~")),data=pheno_geno,family="binomial"))
      
      lm_MAT[which(snps==k),index_1[j]:(index_1[j]+3)]=c(lm_MAT_tmp$coefficients[MAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_MAT_tmp$df.residual)
      lm_PAT[which(snps==k),index_1[j]:(index_1[j]+3)]=c(lm_PAT_tmp$coefficients[PAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_PAT_tmp$df.residual)
      lm_OFF[which(snps==k),index_1[j]:(index_1[j]+3)]=c(lm_OFF_tmp$coefficients[OFF,c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
      
      lm_ALL_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                              sep="~")),data=pheno_geno,family="binomial"))
      
      lm_ALL[which(snps==k),index_2[j]:(index_2[j]+9)]=c(lm_ALL_tmp$coefficients[c(OFF,MAT,PAT),c("Estimate","Std. Error","Pr(>|z|)")],lm_ALL_tmp$df.residual)
    }
    else 
    {
      lm_MAT_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(MAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                             sep="~")),data=pheno_geno))
      lm_PAT_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                             sep="~")),data=pheno_geno))
      lm_OFF_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(OFF,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                             sep="~")),data=pheno_geno))
      
      lm_MAT[which(snps==k),index_1[j]:(index_1[j]+3)]=c(lm_MAT_tmp$coefficients[MAT,c("Estimate","Std. Error","Pr(>|t|)")],lm_MAT_tmp$fstatistic["dendf"])
      lm_PAT[which(snps==k),index_1[j]:(index_1[j]+3)]=c(lm_PAT_tmp$coefficients[PAT,c("Estimate","Std. Error","Pr(>|t|)")],lm_PAT_tmp$fstatistic["dendf"])
      lm_OFF[which(snps==k),index_1[j]:(index_1[j]+3)]=c(lm_OFF_tmp$coefficients[OFF,c("Estimate","Std. Error","Pr(>|t|)")],lm_OFF_tmp$fstatistic["dendf"])
      
      lm_ALL_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                             sep="~")),data=pheno_geno))
      
      lm_ALL[which(snps==k),index_2[j]:(index_2[j]+9)]=c(lm_ALL_tmp$coefficients[c(OFF,MAT,PAT),c("Estimate","Std. Error","Pr(>|t|)")],lm_ALL_tmp$fstatistic["dendf"])
    }
    print(j)
  }
  
}

write.csv(lm_ALL, "lm_r_r_ALL.csv")
write.csv(lm_MAT, "lm_r_r_MAT.csv")
write.csv(lm_OFF, "lm_r_r_PAT.csv")
write.csv(lm_PAT, "lm_r_r_OFF.csv")



