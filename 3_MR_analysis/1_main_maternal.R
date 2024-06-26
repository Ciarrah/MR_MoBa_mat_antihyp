# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(c("BiocGenerics", "Biostrings", "GenomeInfoDb", "GenomicRanges", "IRanges", "Rsamtools", "S4Vectors", "SummarizedExperiment", "VariantAnnotation"))
# remotes::install_github("mrcieu/gwasvcf")
# install.packages("gwasglue-master", repos= NULL, type="source")

# MR and plot

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

# 2. Read in snp-outcome data, append eaf, condition, subclass, substance and restrict to COI.
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

# 3. Append the gene and drug class, check mechanism of action on DrugBank for possible combined instruments per drug subclass.
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

# 4. Check mechanism of action on DrugBank. Analyse Vasodilator antihypertensive drugs separately as they have opposing mechanisms.

full_outcome_2[,"Drug subclass"]=full_outcome_2$drug_subclass
full_outcome_2[which(full_outcome_2$gene=="KCNJ11"),"Drug subclass"]=paste0(full_outcome_2[which(full_outcome_2$gene=="KCNJ11"),"drug_subclass"], 
                                                                            " targeting KCNJ11")
full_outcome_2[which(full_outcome_2$gene=="EDNRA"),"Drug subclass"]=paste0(full_outcome_2[which(full_outcome_2$gene=="EDNRA"),"drug_subclass"], 
                                                                           " targeting EDNRA")
hypertension_exposure=merge(hypertension_exposure, full_outcome_2[,c("SNP", "Drug subclass")], by="SNP")

# 5. Run MR analyses.

results_DC=results_wald_DC=NULL
results_DC=data.frame()
results_wald_DC=data.frame()

for (j in outcomes_of_interest) {
  snp_outcome=format_data(dat = full_outcome_2, type = "outcome", snps = hypertension_exposure$SNP, beta_col = j,
                          se_col = gsub("Beta", "SE", j), eaf_col = "Maternal.effect.allele",
                          effect_allele_col = "Effect.allele", pval_col = gsub("Beta", "Pval", j),
                          other_allele_col = "Other.allele", phenotype_col = "Drug subclass")
  
  # Get unique subclasses
  unique_subclasses=unique(hypertension_exposure$`Drug subclass`)
  
  # Iterate over unique subclasses
  for (k in unique_subclasses) {
    temp=subset(hypertension_exposure, `Drug subclass` == k)
    temp_outcome=subset(snp_outcome, outcome == k)
    harm_temp=harmonise_data(temp, temp_outcome)
    harm_temp=harm_temp[!duplicated(harm_temp), ]
    print(c(j, k, nrow(harm_temp)))
    
    if (nrow(harm_temp) == 1) {
      tmp=mr_singlesnp(harm_temp, parameters = default_parameters(), single_method = "mr_wald_ratio")
      tmp=tmp[-grep("All", tmp$SNP), ]
      tmp[, c("b", "se", "p")]=apply(tmp[, c("b", "se", "p")], 2, as.numeric)
      tmp$`Upper CI`=tmp$b + (1.96 * tmp$se)
      tmp$`Lower CI`=tmp$b - (1.96 * tmp$se)
      tmp=cbind(tmp, k)
      colnames(tmp)[which(colnames(tmp) == "k")]="Drug subclass"
      tmp$Condition="Hypertension"
      tmp$Outcome=j
      tmp$snps=paste0(harm_temp$SNP, collapse = " ")
      results_wald_DC=rbind(results_wald_DC, tmp)
    } else {
      tmp=mr(harm_temp)
      tmp[, c("b", "se", "pval")]=apply(tmp[, c("b", "se", "pval")], 2, as.numeric)
      tmp$`Upper CI`=tmp$b + (1.96 * tmp$se)
      tmp$`Lower CI`=tmp$b - (1.96 * tmp$se)
      tmp=cbind(tmp, k)
      colnames(tmp)[which(colnames(tmp) == "k")]="Drug subclass"
      tmp$Condition="Hypertension"
      tmp$Outcome=j
      tmp$snps=paste0(harm_temp$SNP, collapse = " ")
      results_DC=rbind(results_DC, tmp)
    }
  }
}

results_wald_DC$Outcome=gsub("_MBRN","",gsub("Beta_Mat_","",results_wald_DC$Outcome))
results_DC$Outcome=gsub("_MBRN","",gsub("Beta_Mat_","",results_DC$Outcome))
results_DC$`Drug subclass`=str_to_sentence(results_DC$`Drug subclass`)

# 6. Tidy variables for plot

results_DC$`Drug subclass`=str_wrap(results_DC$`Drug subclass`, width=20)
outcomes_of_interest=gsub("_MBRN","",gsub("Beta_Mat_","",outcomes_of_interest))
results_DC$Outcome[grep("bwt", results_DC$Outcome)]="Birthweight (100g)"
results_DC$Outcome[grep("GA_dys", results_DC$Outcome)]="Gestational age (days)"
results_DC$Outcome[grep("GA_wks", results_DC$Outcome)]="Gestational age (weeks)"
results_DC$Outcome[grep("any_hyp_prg", results_DC$Outcome)]="Hypertensive disorders of pregnancy"
results_DC$Outcome[grep("brth_lngth", results_DC$Outcome)]="Birth length (cm)"
results_DC$Outcome[grep("head_circum", results_DC$Outcome)]="Head circumference (cm)"
results_DC$Outcome[grep("apgar_1min", results_DC$Outcome)]="Apgar score, 1 minute"
results_DC$Outcome[grep("apgar_5min", results_DC$Outcome)]="Apgar score, 5 minutes"
results_DC$Outcome[grep("any_development_delay", results_DC$Outcome)]="Developmental delay"
results_DC$Outcome[grep("prorated_scores", results_DC$Outcome)]="Developmental score"
colnames(results_DC)[which(colnames(results_DC)=="b")]=c("Beta coefficient")

results_wald_DC$`Drug subclass`=str_wrap(results_wald_DC$`Drug subclass`, width=20)
results_wald_DC$Outcome[grep("bwt", results_wald_DC$Outcome)]="Birthweight (100g)"
results_wald_DC$Outcome[grep("GA_dys", results_wald_DC$Outcome)]="Gestational age (days)"
results_wald_DC$Outcome[grep("GA_wks", results_wald_DC$Outcome)]="Gestational age (weeks)"
results_wald_DC$Outcome[grep("any_hyp_prg", results_wald_DC$Outcome)]="Hypertensive disorders of pregnancy"
results_wald_DC$Outcome[grep("brth_lngth", results_wald_DC$Outcome)]="Birth length (cm)"
results_wald_DC$Outcome[grep("head_circum", results_wald_DC$Outcome)]="Head circumference (cm)"
results_wald_DC$Outcome[grep("apgar_1min", results_wald_DC$Outcome)]="Apgar score, 1 minute"
results_wald_DC$Outcome[grep("apgar_5min", results_wald_DC$Outcome)]="Apgar score, 5 minutes"
results_wald_DC$Outcome[grep("any_development_delay", results_wald_DC$Outcome)]="Developmental delay"
results_wald_DC$Outcome[grep("prorated_scores", results_wald_DC$Outcome)]="Developmental score"
colnames(results_wald_DC)[which(colnames(results_wald_DC)=="b")]=c("Beta coefficient")
results_wald_DC$`Drug subclass`=str_to_sentence(results_wald_DC$`Drug subclass`)

colnames(results_wald_DC)[which(colnames(results_wald_DC)=="p")]="pval"

results_ivw_DC=results_DC[which(results_DC$method=="Inverse variance weighted"),]

results_ivw_DC$outcome=str_to_sentence(results_ivw_DC$outcome)
results_ivw_DC$`MR method`="IVW"
results_wald_DC$`MR method`="Wald ratio"
results_wald_DC$nsnp=1

saved_wald_results=results_wald_DC
saved_ivw_results=results_ivw_DC

results_all_DC=rbind(saved_ivw_results[,c("Beta coefficient","se","pval","Upper CI","Lower CI","Drug subclass","MR method","Outcome","snps","nsnp")],
                     saved_wald_results[which(saved_wald_results$`Drug subclass`%!in%saved_ivw_results$outcome),c("Beta coefficient","se","pval","Upper CI","Lower CI","Drug subclass","MR method","Outcome","snps","nsnp")])

results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="Wald ratio")]=paste0(results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="Wald ratio")]," (1 SNP)")
results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="IVW")]=paste0(results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="IVW")]," (",
                                                                                results_all_DC$nsnp[which(results_all_DC$`MR method`=="IVW")], " SNPs)")
results_all_DC$Condition="Hypertension"

# Rescale birthweight to per 100g increase

results_all_DC[which(results_all_DC$Outcome=="Birthweight (100g)"),c("Beta coefficient","Upper CI", "Lower CI")]=
  results_all_DC[which(results_all_DC$Outcome=="Birthweight (100g)"),c("Beta coefficient","Upper CI", "Lower CI")]/100

# Ensure numeric. Exponentiate coefficients from binary variables

results_all_DC[,c("Beta coefficient", "se", "pval","Upper CI", "Lower CI")]=apply(results_all_DC[,c("Beta coefficient", "se", "pval","Upper CI", "Lower CI")], 2, as.numeric)
results_all_DC[which(results_all_DC$Outcome=="Hypertensive disorders of pregnancy"),c("Beta coefficient", "Upper CI", "Lower CI")]=
  exp(results_all_DC[which(results_all_DC$Outcome=="Hypertensive disorders of pregnancy"),c("Beta coefficient", "Upper CI", "Lower CI")])

# 7. Transform variables from per sd increase to 10mmHg for interpretability. Flip the sign of the MR estimates so 
# they are representative of drug effects (decreasing). Divide by SD and multiply by 10 to interpret on per 10mmHg scale.

SD=19.38 # Biobank showcase value

results_all_DC$`Beta coefficient`=(results_all_DC$`Beta coefficient`/19.38)*10
results_all_DC$`Beta coefficient`=-results_all_DC$`Beta coefficient`

# transform and swap CI
results_all_DC$`Upper CI`=(results_all_DC$`Upper CI`/19.38)*10

results_all_DC$`Lower CI`=(results_all_DC$`Lower CI`/19.38)*10

tmp=results_all_DC$`Upper CI`
results_all_DC$`Upper CI`=results_all_DC$`Lower CI`
results_all_DC$`Lower CI`=tmp


# 8. Plot results, "Preeclampsia" panel with modified x-axis limits.

results_all_DC$Outcome=str_wrap(results_all_DC$Outcome, width=20)
write.csv(results_all_DC ,"full_maternal_results.csv")

results_all_DC=read.csv("full_maternal_results.csv")
colnames(results_all_DC)=gsub("\\.", " ", colnames(results_all_DC))

filename=paste0("hypertension_subclass_maternal_r&r.tiff")
plot_normal_scale=ggplot(results_all_DC[results_all_DC$Outcome != "Hypertensive\ndisorders of\npregnancy", ], aes(x = `Beta coefficient`, y = `Drug subclass`)) +
  geom_point(size = 2, na.rm = TRUE) +
  scale_y_discrete(position = "left") +
  geom_linerange(aes(xmin = `Lower CI`, xmax = `Upper CI`), na.rm = TRUE) +
  facet_grid2(`Condition` ~ `Outcome`, scales = "free", space = "free") +
  force_panelsizes(cols = c(rep(0.5, 7))) +
  geom_vline(data = subset(results_all_DC, Outcome=="Hypertensive\ndisorders of\npregnancy"), aes(xintercept = 0)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        strip.text.y = element_blank(),
        axis.text.y.right = element_blank(),
        axis.text.y=element_blank())

plot_gest_hypertension=ggplot(results_all_DC[results_all_DC$Outcome%in%c("Hypertensive\ndisorders of\npregnancy"), ], aes(x = `Beta coefficient`, y = `Drug subclass`)) +
  geom_point(size = 2) +
  scale_y_discrete(position = "left") +
  scale_x_continuous(trans = "log10")+
  geom_linerange(aes(xmin = `Lower CI`, xmax = `Upper CI`)) +
  facet_grid2(`Condition` ~ `Outcome`, scales = "free", space = "free")+
  geom_vline(data = subset(results_all_DC, Outcome%in% c("Hypertensive\ndisorders of\npregnancy")), aes(xintercept = 1)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines"),
        axis.text.y.left = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 15),
        axis.title.x = element_blank(),
        text = element_text(size = 12),
        strip.text.y = element_blank(),
        axis.text.y.right = element_blank())+
  guides(color = "none")
x_title="Odds Ratio (OR), binary outcome, or beta, continuous outcome, and 95% confidence interval (CI) for a 10mmHg decrease in systolic blood pressure"
tiff(filename, units="cm", width=58, height=22, res=300)
grid.arrange(plot_gest_hypertension, plot_normal_scale,ncol = 2, widths = c(1, 4), bottom=x_title)
dev.off()

# Results files

results_table_Maternal_DC=as.data.frame(matrix(0,ncol=10, nrow=(nrow(results_all_DC))))
colnames(results_table_Maternal_DC)=c("Outcome","Beta coefficient (95% CI)", "Standard error", "Pvalue", "No. SNPs", "Drug subclass", "Condition", "MR method", "Gene", "rsID")
results_table_Maternal_DC[,"Outcome"]=results_all_DC$Outcome
results_table_Maternal_DC[,"Beta coefficient (95% CI)"]=c(paste0(round(results_all_DC$`Beta coefficient`, 2), " (", round(results_all_DC$`Lower CI`,2), ", ", round(results_all_DC$`Upper CI`,2), ")"))
results_table_Maternal_DC[,"Standard error"]=round(results_all_DC$se,2)
results_table_Maternal_DC[,"MR method"]=results_all_DC$`MR method`
results_table_Maternal_DC[,"Pvalue"]=round(results_all_DC$pval,2)
results_table_Maternal_DC[,"Drug subclass"]=results_all_DC$`Drug subclass`
results_table_Maternal_DC[,"Drug subclass"]=gsub("\n"," ",results_all_DC$`Drug subclass`)
results_table_Maternal_DC[,"No. SNPs"]=results_all_DC$nsnp
results_table_Maternal_DC[,"Gene"]=results_all_DC$Gene
results_table_Maternal_DC[,"rsID"]=results_all_DC$snps
results_table_Maternal_DC=results_table_Maternal_DC[-which(results_table_Maternal_DC$Outcome=="Developmental delay"),]
write.csv(results_table_Maternal_DC,"hypertension_subclass_maternal.csv" )

write.csv(results_all_DC, "gest_r_full_file_mat.csv")
