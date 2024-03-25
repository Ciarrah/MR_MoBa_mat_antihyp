# Parallise over snps, run on cluster

library(data.table)
`%!in%` = Negate(`%in%`)

args=commandArgs(trailingOnly=T)
snp=(args[1])
print(snp)
#snp=gsub("(^\\\"|\\\"$)", "", snp)
snp=snp[1]
#snp=strsplit(gsub("[^[:alnum:] ]", "", snp), " +")[[1]]
paste0(c("TMP",snp))

# Load in the phenotypic and genotypic dataset. Load in the list of SNPs for each drug and condition.
# Column names are the first SNP, rbind the column name for each condition to include the SNP.

pheno_geno=fread("pheno_geno_hyp.txt")

# 1. Phenotype=offspring_snp+covariates
# 2. Phenotype=father_snp+covariates
# 3. Phenotype=mother_snp+covariates
# 4. Phenotype=offspring_snp+father_snp+mother_snp+covariates
# Run SNP-wise regression of each SNP on each outcome of interest.

# Then run a regression of all the relevant SNPs for each drug subclass on the outcome of interest
# Include PCs and batch variables I.e. one SNP per regression. And then run regression separately across all the SNPs.                                               
                                                    
outcome_list=grep("MBRN",colnames(pheno_geno), value=T)[c(3,4,6,10,12,13,14,15,16,26,27,28,32,33,34,35,36,37)]

# Load in SNP details, create an output table of coefficients, SE, pval, DF for each outcome
# with eaf, other, and allele freq appended.

# less SNPs are available in SNP_details than super_set because they may not be available in MoBa

lm_MAT=lm_PAT=lm_OFF=as.data.frame(matrix(NA,nrow=1,ncol=(4*length(outcome_list))))
lm_ALL=as.data.frame(matrix(NA,nrow=1,ncol=(10*length(outcome_list))))
MAT=PAT=OFF=ALL=NULL
rownames(lm_MAT)=rownames(lm_PAT)=rownames(lm_OFF)=rownames(lm_ALL)=paste0(snp)
index_1=seq(from=1,72,by=4)
colnames(lm_MAT)=colnames(lm_PAT)=colnames(lm_OFF)=paste0(rep(c("Beta","SE","Pval","DF"),length(outcome_list)),"_",rep(outcome_list,each=4))
index_2=seq(from=1,180,by=10)
colnames(lm_ALL)=paste0(c(paste0(c(rep("Beta",3),rep("SE",3),rep("Pval",3)),rep(c("_Off","_Mat","_Pat"),3)),"DF"),"_",rep(outcome_list,each=10))

# lm for continuous models or categorical ones
# glm for binary

# Ensure correct classes

pheno_geno$preeclmps_MBRN[which(pheno_geno$preeclmps_MBRN%in%c(1:3))]=1
pheno_geno$preeclmps_MBRN=as.factor(pheno_geno$preeclmps_MBRN)
pheno_geno$hyptnsn_pre_prg_MBRN=as.factor(pheno_geno$hyptnsn_pre_prg_MBRN)
pheno_geno$hyptnsn_prg_MBRN=as.factor(pheno_geno$hyptnsn_prg_MBRN)
pheno_geno$eclmps_MBRN=as.factor(pheno_geno$eclmps_MBRN)
pheno_geno$ealry_preeclmps_MBRN=as.factor(pheno_geno$early_preeclmps_MBRN)
pheno_geno$early_eclmps_MBRN=as.factor(pheno_geno$early_eclmps_MBRN)
pheno_geno$HELLP_MBRN=as.factor(pheno_geno$HELLP_MBRN)
pheno_geno$dlvry_initin_MBRN=as.factor(pheno_geno$dlvry_initin_MBRN)
pheno_geno$dlvry_initin_MBRN[which(pheno_geno$dlvry_initin_MBRN=="NA")]=NA

# Create a loop that runs linear regression for each continuous outcome, logistic for binary variables and 
# continuous for categorical variables with more than one variable.

# Loop in bash script over each SNP, pass in SNP to R and output each file of all regression results to csv.

snps=grep(paste0(snp),colnames(pheno_geno),value=T)
MAT=grep("mothers",snps,value=T)
PAT=grep("fathers",snps,value=T)
OFF=grep("offspring",snps,value=T)

for (j in 1:length(outcome_list))
{
  if ((is.factor(pheno_geno[[outcome_list[j]]]))==T&&length(levels(pheno_geno[[outcome_list[j]]]))>2)
  {
    lm_MAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(MAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_PAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_OFF_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_ALL_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    
    lm_MAT[1,index_1[j]:(index_1[j]+3)]=c(lm_MAT_tmp$coefficients[MAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_MAT_tmp$df.residual)
    lm_PAT[1,index_1[j]:(index_1[j]+3)]=c(lm_PAT_tmp$coefficients[PAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_PAT_tmp$df.residual)
    lm_OFF[1,index_1[j]:(index_1[j]+3)]=c(lm_OFF_tmp$coefficients[OFF,c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
    
    lm_ALL_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    
    lm_ALL[1,index_2[j]:(index_2[j]+9)]=c(lm_ALL_tmp$coefficients[c(OFF,MAT,PAT),c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
  }
  else if ((is.factor(pheno_geno[[outcome_list[j]]]))==T&&length(levels(pheno_geno[[outcome_list[j]]]))==2)
  {
    lm_MAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(MAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_PAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_OFF_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    
    lm_MAT[1,index_1[j]:(index_1[j]+3)]=c(lm_MAT_tmp$coefficients[MAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_MAT_tmp$df.residual)
    lm_PAT[1,index_1[j]:(index_1[j]+3)]=c(lm_PAT_tmp$coefficients[PAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_PAT_tmp$df.residual)
    lm_OFF[1,index_1[j]:(index_1[j]+3)]=c(lm_OFF_tmp$coefficients[OFF,c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
    
    lm_ALL_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    
    lm_ALL[1,index_2[j]:(index_2[j]+9)]=c(lm_ALL_tmp$coefficients[c(OFF,MAT,PAT),c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
  }
  
  else if ((is.factor(pheno_geno[[outcome_list[j]]]))==T&&length(levels(pheno_geno[[outcome_list[j]]]))==1)
  {
    lm_MAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(MAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_PAT_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    lm_OFF_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    
    lm_MAT[1,index_1[j]:(index_1[j]+3)]=c(lm_MAT_tmp$coefficients[MAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_MAT_tmp$df.residual)
    lm_PAT[1,index_1[j]:(index_1[j]+3)]=c(lm_PAT_tmp$coefficients[PAT,c("Estimate","Std. Error","Pr(>|z|)")],lm_PAT_tmp$df.residual)
    lm_OFF[1,index_1[j]:(index_1[j]+3)]=c(lm_OFF_tmp$coefficients[OFF,c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
    
    lm_ALL_tmp=summary(glm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                            sep="~")),data=pheno_geno,family="binomial"))
    
    lm_ALL[1,index_2[j]:(index_2[j]+9)]=c(lm_ALL_tmp$coefficients[c(OFF,MAT,PAT),c("Estimate","Std. Error","Pr(>|z|)")],lm_OFF_tmp$df.residual)
    
  }
  else 
  {
    lm_MAT_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(MAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                           sep="~")),data=pheno_geno))
    lm_PAT_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                           sep="~")),data=pheno_geno))
    lm_OFF_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(OFF,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                           sep="~")),data=pheno_geno))
    
    lm_MAT[1,index_1[j]:(index_1[j]+3)]=c(lm_MAT_tmp$coefficients[MAT,c("Estimate","Std. Error","Pr(>|t|)")],lm_MAT_tmp$fstatistic["dendf"])
    lm_PAT[1,index_1[j]:(index_1[j]+3)]=c(lm_PAT_tmp$coefficients[PAT,c("Estimate","Std. Error","Pr(>|t|)")],lm_PAT_tmp$fstatistic["dendf"])
    lm_OFF[1,index_1[j]:(index_1[j]+3)]=c(lm_OFF_tmp$coefficients[OFF,c("Estimate","Std. Error","Pr(>|t|)")],lm_OFF_tmp$fstatistic["dendf"])
    
    lm_ALL_tmp=summary(lm(as.formula(paste(outcome_list[j],paste(c(OFF,MAT,PAT,grep("age|sx|PC|batch",colnames(pheno_geno),value=TRUE)),collapse="+"),
                                           sep="~")),data=pheno_geno))
    
    lm_ALL[1,index_2[j]:(index_2[j]+9)]=c(lm_ALL_tmp$coefficients[c(OFF,MAT,PAT),c("Estimate","Std. Error","Pr(>|t|)")],lm_ALL_tmp$fstatistic["dendf"])
  }
  print(j)
}

## Export out all files and append index to filename.

head(lm_MAT)
write.csv(lm_MAT,paste0("hyp_lm_MAT_",paste0(snp),".csv"))
write.csv(lm_PAT,paste0("hyp_lm_PAT_",paste0(snp),".csv"))
write.csv(lm_OFF,paste0("hyp_lm_OFF_",paste0(snp),".csv"))
write.csv(lm_ALL,paste0("hyp_lm_ALL_",paste0(snp),".csv"))
