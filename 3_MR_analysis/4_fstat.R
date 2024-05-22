library(readr);library(TwoSampleMR);library(data.table);library(tidyr);library(stringr);library(ggplot2);library(ggh4x);library(vcfR)
library(gwasvcf);library(gwasglue);library(ggplot2);library(gridExtra);library(readxl);library(LDlinkR);library(readxl);library(dplyr)

`%!in%` = Negate(`%in%`)
options(scipen = 999)

## 1. Check the overlap between the snp-exposure GWAS and snp-outcome GWAS

# Load full, unclumped SNP-exposure gwas with all significance level
hypertension_full=gwasvcf::query_gwas("ukb-b-20175.vcf.gz", pval = 1)
hypertension_full=gwasglue::gwasvcf_to_TwoSampleMR(hypertension_full, "exposure")

# Append gene info from gencode, subset to SNPs with RSID
b37=fread("common_all_20180423.vcf.gz")
tmp=hypertension_full[hypertension_full$SNP%in%b37$ID,]
tmp=merge(tmp,b37[,c("ID","INFO")], by.x="SNP",by.y="ID")
tmp$gene=gsub('.*GENEINFO=([^;:]+).*', '\\1', tmp[["INFO"]])

# List of snps in MoBa
moba_snps=read_excel("MoBa_SNP_list.xlsx", col_names = FALSE)
colnames(moba_snps)[1]="RSID"
moba_snps$RSID=gsub("_.","",moba_snps$RSID)
hypertension_exposure=tmp[which(tmp$SNP%in%moba_snps$RSID),]
hypertension_exposure=hypertension_exposure[which(hypertension_exposure$pval.exposure<5*10^-8),]
hypertension_exposure=clump_data(hypertension_exposure, clump_r2=0.01)
# 23,432 SNPs with GENCODE info, 292 SNPs at GWAS significance level, 7 remain after clumping

## 2. Read in snp-outcome data, append eaf, condition, subclass, substance and restrict to COI.

outcome=read_excel("Regression_results.xlsx")
colnames(outcome)[1]="RSID"
eaf=read.csv("summary_table_genetic_trio.csv")[,-1]
full_outcome=merge(outcome, eaf, by="RSID")
colnames(full_outcome)[1]="SNP"
full_outcome$CHR=as.character(full_outcome$CHR)
rm(outcome, eaf)

# Check no SNPs missing rsID
no_rsid=dim(full_outcome[which(is.na(full_outcome$CHR)),"SNP"])
print(no_rsid)

## 3. Append the gene and drug class, check mechanism of action on DrugBank for possible combined instruments per drug subclass.

full_outcome$gene=NULL
double_check=as.data.frame(matrix(NA,nrow=0, ncol=4))
gene_pos=read_delim("gene_chr_pos_hyp.txt", col_names = F)
colnames(gene_pos)=c("chr.exposure","pos.exposure","pos.end.exposure","gene")
for (i in 1:nrow(full_outcome))
{
  tmp_2=gene_pos[which(gene_pos$chr.exposure==full_outcome[i,"CHR"]),]
  
  # For chromosome match, is the bp in the range
  gene=tmp_2$gene[which(tmp_2$pos.exposure <= full_outcome[i,"POS"] & full_outcome[i,"POS"] <= tmp_2$pos.end.exposure)]
  if (length(gene)==1)
  {
    full_outcome[i,"gene"]=tmp_2$gene[which(tmp_2$pos.exposure <= full_outcome[i,"POS"] & full_outcome[i,"POS"] <= tmp_2$pos.end.exposure)]
  }
  else if (length(gene)>1) # check no overlap -- snp assigned two genes
  {
    double_check=rbind(double_check, full_outcome[i,grep("SNP|CHR|gene|POS", colnames(full_outcome))])
    print(tmp_2$gene[which(tmp_2$pos.exposure <= full_outcome[i,"POS"] & full_outcome[i,"POS"] <= tmp_2$pos.end.exposure)])
  }
}

genetic_targets=read.csv("identified_genes_hypertension_190523.csv")
hypertension_bnf_code_to_drug_class=read_csv("hypertension_bnf_code_to_drug_class.csv")

genetic_targets=merge(genetic_targets, hypertension_bnf_code_to_drug_class, all=T)
genetic_targets[grep("2040", genetic_targets$bnf_code),"drug_subclass"]="Beta-adrenoceptor blocking drugs"
genetic_targets[grep("20501", genetic_targets$bnf_code),"drug_subclass"]="Vasodilator antihypertensive drugs"
genetic_targets[grep("20502", genetic_targets$bnf_code),"drug_subclass"]="Centrally-acting antihypertensive drugs"
genetic_targets[grep("20503", genetic_targets$bnf_code),"drug_subclass"]="Adrenergic neurone blocking drugs"
genetic_targets[grep("20504", genetic_targets$bnf_code),"drug_subclass"]="Alpha-adrenoceptor blocking drugs"
genetic_targets[grep("20505", genetic_targets$bnf_code),"drug_subclass"]="Renin-angiotensin system drugs"
genetic_targets=genetic_targets[!is.na(genetic_targets$substance),]
genetic_targets=distinct(genetic_targets[which(genetic_targets$gene%in%full_outcome$gene),colnames(genetic_targets)!="substance"])
full_outcome_2=distinct(merge(full_outcome, genetic_targets[,c("gene", "drug_subclass")], by="gene"))

# Check all SNPs assigned a BNF code, drug subclass and gene
dim(full_outcome_2$bnf_code[is.na(full_outcome_2$bnf_code)])
dim(full_outcome_2$gene[is.na(full_outcome_2$gene)])
dim(full_outcome_2$drug_subclass[is.na(full_outcome_2$drug_subclass)])

outcomes_of_interest=colnames(full_outcome_2)[grep("Beta_Mat", colnames(full_outcome_2))]
outcomes_of_interest=outcomes_of_interest[grep("apgar|GA_dys|bwt|lngth|circum|any_hyp|any_develop|prorated_score",
                                               outcomes_of_interest)]

## 4. Check mechanism of action on DrugBank. Analyse Vasodilator antihypertensive drugs separately as they have opposing mechanisms.

full_outcome_2[,"Drug subclass"]=full_outcome_2$drug_subclass
full_outcome_2[which(full_outcome_2$gene=="KCNJ11"),"Drug subclass"]=paste0(full_outcome_2[which(full_outcome_2$gene=="KCNJ11"),"drug_subclass"], 
                                                                            " targeting KCNJ11")
full_outcome_2[which(full_outcome_2$gene=="EDNRA"),"Drug subclass"]=paste0(full_outcome_2[which(full_outcome_2$gene=="EDNRA"),"drug_subclass"], 
                                                                           " targeting EDNRA")
hypertension_exposure=merge(hypertension_exposure, full_outcome_2[,c("SNP", "Drug subclass")], by="SNP")

# 5. Fstatistic

# Calculate pre-fstat
hypertension_exposure$pre_fstat=(hypertension_exposure$beta.exposure^2)/(hypertension_exposure$se.exposure^2)

# Take mean across substance/subclass
df1=hypertension_exposure[,c("pre_fstat", "Drug subclass")]
df1=df1[!duplicated(df1),]
fstat_subclass=tapply(df1$pre_fstat, df1$`Drug subclass`, mean)

write.csv(as.data.frame(fstat_subclass), "hypertension_fstat_subclass.csv",
          row.names=T)
