## R&R MR analysis

library(readr);library(TwoSampleMR);library(data.table);library(tidyr);library(stringr);library(ggplot2);library(ggh4x);library(vcfR);library(grid)
library(gwasvcf);library(gwasglue);library(gridExtra);library(readxl);library(LDlinkR);library(readxl);library(dplyr);library(openxlsx);library(openxlsx);library(gridExtra)
`%!in%` = Negate(`%in%`)
options(scipen = 999)

## 1. Check the overlap between the snp-exposure GWAS and snp-outcome GWAS

# Load full, unclumped SNP-exposure gwas with all significance level
hypertension_full=gwasvcf::query_gwas("ukb-b-20175.vcf.gz", pval = 1)
hypertension_full=gwasglue::gwasvcf_to_TwoSampleMR(hypertension_full, "exposure") #9797146 without dupliate snps

# Append gene info from gencode, subset to SNPs with RSID
b37=fread("common_all_20180423.vcf.gz")
tmp=hypertension_full[hypertension_full$SNP%in%b37$ID,] #9291736
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


# revised outcomes

outcome_original=read_excel("Regression_results.xlsx")
outcome_additional=read.csv("lm_r_r_ALL.csv")
colnames(outcome_original)[1]=colnames(outcome_additional)[1]="SNP"
outcome=as.data.frame(outcome_additional)
outcome_original=as.data.frame(outcome_original[,grep("SNP|Pat_GA_wks_MBRN|Pat_brth_lngth_MBRN|any_hyp|Pat_head_circum_MBRN|Pat_apgar_1min_MBRN|Pat_apgar_5min_MBRN|Pat_bwt_MBRN",
                                                      colnames(outcome_original))])
outcome=merge(outcome, outcome_original, by="SNP")

colnames(outcome)[1]="RSID"
eaf=read.csv("summary_table_genetic_trio_hyp.csv")[,-1]
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

genetic_targets=read.csv("identified_genes_hypertension.csv")
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

outcomes_of_interest=colnames(full_outcome_2)[grep("Beta_Pat", colnames(full_outcome_2))]

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
                          se_col = gsub("Beta", "SE", j), eaf_col = "Paternal.effect.allele",
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

results_wald_DC$Outcome=gsub("_MBRN","",gsub("Beta_Pat_","",results_wald_DC$Outcome))
results_DC$Outcome=gsub("_MBRN","",gsub("Beta_Pat_","",results_DC$Outcome))
results_DC$`Drug subclass`=str_to_sentence(results_DC$`Drug subclass`)

# 6. Tidy variables for plot

results_DC$`Drug subclass`=str_wrap(results_DC$`Drug subclass`, width=20)
results_DC$Outcome[grep("GA_dys", results_DC$Outcome)]="Gestational age (days)"
results_DC$Outcome[grep("any_hyp_prg", results_DC$Outcome)]="Hypertensive disorders of pregnancy"
results_DC$Outcome[grep("brth_lngth", results_DC$Outcome)]="Birth length (cm)"
results_DC$Outcome[grep("head_circum", results_DC$Outcome)]="Head circumference (cm)"
results_DC$Outcome[grep("apgar_1min", results_DC$Outcome)]="Apgar score, 1 minute"
results_DC$Outcome[grep("apgar_5min", results_DC$Outcome)]="Apgar score, 5 minutes"
results_DC$Outcome[grep("prorated_scores", results_DC$Outcome)]="Developmental score"
results_DC$Outcome[grep("BW_GA", results_DC$Outcome)]="Birthweight z-score"
results_DC$Outcome[grep("Congenital_malf", results_DC$Outcome)]="Congenital malformation"

colnames(results_DC)[which(colnames(results_DC)=="b")]=c("Beta coefficient")

results_wald_DC$`Drug subclass`=str_wrap(results_wald_DC$`Drug subclass`, width=20)
results_wald_DC$Outcome[grep("GA_dys", results_wald_DC$Outcome)]="Gestational age (days)"
results_wald_DC$Outcome[grep("any_hyp_prg", results_wald_DC$Outcome)]="Hypertensive disorders of pregnancy"
results_wald_DC$Outcome[grep("brth_lngth", results_wald_DC$Outcome)]="Birth length (cm)"
results_wald_DC$Outcome[grep("head_circum", results_wald_DC$Outcome)]="Head circumference (cm)"
results_wald_DC$Outcome[grep("apgar_1min", results_wald_DC$Outcome)]="Apgar score, 1 minute"
results_wald_DC$Outcome[grep("apgar_5min", results_wald_DC$Outcome)]="Apgar score, 5 minutes"
results_wald_DC$Outcome[grep("prorated_scores", results_wald_DC$Outcome)]="Developmental score"
results_wald_DC$Outcome[grep("BW_GA", results_wald_DC$Outcome)]="Birthweight z-score"
results_wald_DC$Outcome[grep("Congenital_malf", results_wald_DC$Outcome)]="Congenital malformation"

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

# 7. Transform variables from per sd increase to 10mmHg for interpretability. Flip the sign of the MR estimates so 
# they are representative of drug effects (decreasing). Divide by SD and multiply by 10 to interpret on per 10mmHg scale.

results_all_DC=rbind(saved_ivw_results[,c("Beta coefficient","se","pval","Upper CI","Lower CI","Drug subclass","MR method","Outcome","snps","nsnp")],
                     saved_wald_results[which(saved_wald_results$`Drug subclass`%!in%saved_ivw_results$outcome),c("Beta coefficient","se","pval","Upper CI","Lower CI","Drug subclass","MR method","Outcome","snps","nsnp")])

results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="Wald ratio")]=paste0(results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="Wald ratio")]," (1 SNP)")
results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="IVW")]=paste0(results_all_DC$`Drug subclass`[which(results_all_DC$`MR method`=="IVW")]," (",
                                                                                results_all_DC$nsnp[which(results_all_DC$`MR method`=="IVW")], " SNPs)")
results_all_DC$Condition="Hypertension"

results_all_DC$`Beta coefficient`=(results_all_DC$`Beta coefficient`/19.38)*10
results_all_DC$`Beta coefficient`=-results_all_DC$`Beta coefficient`

# transform and swap CI
results_all_DC$`Upper CI`=-(results_all_DC$`Upper CI`/19.38)*10

results_all_DC$`Lower CI`=-(results_all_DC$`Lower CI`/19.38)*10

results_all_DC$tmp=results_all_DC$`Upper CI`
results_all_DC$`Upper CI`=results_all_DC$`Lower CI`
results_all_DC$`Lower CI`=results_all_DC$tmp

results_all_DC[which(results_all_DC$Outcome%in%c("Congenital malformation", "Hypertensive disorders of pregnancy")),c("Beta coefficient", "Upper CI", "Lower CI")]=
  exp(results_all_DC[which(results_all_DC$Outcome%in%c("Congenital malformation", "Hypertensive disorders of pregnancy")),c("Beta coefficient", "Upper CI", "Lower CI")])

results_all_DC$Outcome=str_wrap(results_all_DC$Outcome, width=20)

# 8. Plot results, "Preeclampsia" panel with modified x-axis limits.

# write.xlsx(results_all_DC ,"full_paternal_results.xlsx")

results_all_DC=read.xlsx("full_paternal_results.xlsx")
colnames(results_all_DC)=gsub("\\.", " ", colnames(results_all_DC))

# Post review changes

results_all_DC$`Drug subclass 2`=NA
results_all_DC$`Drug subclass 2`[grep("Calcium-channel", results_all_DC$`Drug subclass`)]="Calcium-channel blockers targeting CACNB2 (3 SNPs)"
results_all_DC$`Drug subclass 2`[grep("Ednra", results_all_DC$`Drug subclass`)]="Vasodilator Antihypertensive Drugs targeting EDNRA (1 SNP)"
results_all_DC$`Drug subclass 2`[grep("Beta-adrenoceptor", results_all_DC$`Drug subclass`)]="Beta-adrenoceptor Blocking drugs targeting ADRB1 (1 SNP)"
results_all_DC$`Drug subclass 2`[grep("Potassium-sparing", results_all_DC$`Drug subclass`)]="Potassium-sparing Diuretics and Aldosterone Antagonists targeting SCNN1D (1 SNP)"
results_all_DC$`Drug subclass 2`[grep("Kcnj11", results_all_DC$`Drug subclass`)]="Vasodilator Antihypertensive Drugs targeting KCNJ11 (1 SNP)"

results_all_DC$tmp=results_all_DC$`Drug subclass`
results_all_DC$`Drug subclass`=results_all_DC$`Drug subclass 2`
results_all_DC$`Drug subclass`=sapply(strwrap(results_all_DC$`Drug subclass`, width = 20, simplify = FALSE), paste, collapse = "\n")

# Format for table

sigfigs <- function(x, n_sigfig=3){
  abs_x = abs(x)
  if (abs_x > 0 & abs_x < 1)
  {
    n_sigfig=n_sigfig - 1
  }
  w=signif(x, n_sigfig)
  if (nchar(w)==4 && substr(w, 1, 1) == "0")
  {
    w=paste0(w, 0)
  }
  if (nchar(w)==3 && substr(w, 1, 1) == "0")
  {
    w=paste0(w, 0, 0)
  }
  if (nchar(w)==4 && substr(w, 1, 2) == "-0")
  {
    w=paste0(w, 0, 0)
  }
  if (nchar(w)==5 && substr(w, 1, 2) == "-0")
  {
    w=paste0(w, 0)
  }
  return(w)
}


results_all_DC$f_beta=lapply(results_all_DC$`Beta coefficient`, sigfigs)
results_all_DC$f_lower=lapply(results_all_DC$`Lower CI`, sigfigs)
results_all_DC$f_upper=lapply(results_all_DC$`Upper CI`, sigfigs)

results_all_DC$'Estimate (95% CI)'=paste0(results_all_DC$f_beta, " (",
                                          gsub(" ", "", results_all_DC$f_lower), ", ",
                                          results_all_DC$f_upper, ")")

results_all_DC$`Drug subclass`=gsub("\n", " ", results_all_DC$`Drug subclass`)

for (k in unique(results_all_DC$Outcome))
{
  
  filtered_results <- results_all_DC[results_all_DC$Outcome %in% c(k), ]
  
  # Table
  
  data_table_1 <- ggplot(filtered_results, aes(y = `Drug subclass`)) +
    geom_text(aes(x = 0, label = `Drug subclass`), hjust = 0) +
    geom_text(aes(x = 4, label = `Estimate (95% CI)`), hjust = 1, vjust = 0.5) +
    geom_hline(yintercept = seq(0.5, nlevels(factor(filtered_results$`Drug subclass`)) - 0.5), 
               color = "black", alpha = 0.5) +
    theme_void() +
    coord_cartesian(xlim = c(0, 4)) +
    annotate("text", x = 0.5, y = 5.5, label = "Drug subclass", fontface = "bold") +
    annotate("text", x = 3.5, y = 5.5, label = "Estimate (95% CI)", fontface = "bold") +
    labs(title = NULL) +
    theme(
      strip.background = element_blank(),
      plot.margin = unit(c(1, 2, 1, 2), "lines"),
      text = element_text(size = 9))
  
  # Plot
  
  if (length(grep("Hypertensive|Congenital", k)) == 0) {
    plot_1 <- ggplot(results_all_DC[results_all_DC$Outcome %in% c(k), ], aes(x = `Beta coefficient`, y = `Drug subclass`)) +
      geom_point(size = 2, na.rm = TRUE) +
      scale_y_discrete(position = "left") +
      geom_linerange(aes(xmin = `Lower CI`, xmax = `Upper CI`), na.rm = TRUE) + 
      geom_vline(aes(xintercept = 0)) +
      theme_bw() +
      theme(text = element_text(size=10),
            strip.background = element_blank(),
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(),
            axis.title.x = element_blank(), 
            axis.title.y = element_blank())
    
  } else if (length(grep("Congenital", k)) == 1){
    plot_1 <- ggplot(results_all_DC[results_all_DC$Outcome %in% c(k), ], aes(x = `Beta coefficient`, y = `Drug subclass`)) +
      geom_point(size = 2) +
      scale_y_discrete(position = "left") +
      scale_x_continuous(trans = "log10") +
      geom_linerange(aes(xmin = `Lower CI`, xmax = `Upper CI`)) +
      geom_vline(data = subset(results_all_DC, Outcome %in% c(k)), aes(xintercept = 1)) +
      theme_bw() +
      theme(text = element_text(size=10),
            strip.background = element_blank(),
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(),
            axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
      guides(color = "none")
  }
  else if (length(grep("Hypertensive", k)) == 1){
    plot_1 <- ggplot(results_all_DC[results_all_DC$Outcome %in% c(k), ], aes(x = `Beta coefficient`, y = `Drug subclass`)) +
      geom_point(size = 2) +
      scale_y_discrete(position = "left") +
      scale_x_continuous(trans = "log10") +
      geom_linerange(aes(xmin = `Lower CI`, xmax = `Upper CI`)) +
      geom_vline(data = subset(results_all_DC, Outcome %in% c(k)), aes(xintercept = 1)) +
      theme_bw() +
      theme(text = element_text(size=10),
            strip.background = element_blank(),
            axis.text.y=element_blank(), 
            axis.ticks.y=element_blank(),
            axis.title.x = element_blank(), 
            axis.title.y = element_blank()) +
      guides(color = "none")
  }
  
  assign(paste0("table_plot_",which(unique(results_all_DC$Outcome)==k)),
         grid.arrange(data_table_1, plot_1, ncol = 2, top = textGrob(gsub("\n", " ", k),
                                                                     gp = gpar(fontface = "bold", fontsize=15))))
}

x_title="Beta, continuous outcome, and 95% confidence interval (CI) for a 10mmHg decrease in systolic blood pressure"
tiff(filename=paste0("a_paternal_w_table_1_r&r.tiff"),
     units="cm", width=45, height=45, res=300)
plot_a=grid.arrange(table_plot_9, table_plot_10, table_plot_7,
                    ncol = 1, bottom=x_title)
plot_a
dev.off()

x_title="Odds ratio (OR), binary outcome, or beta, continuous outcome, and 95% confidence interval (CI) for a 10mmHg decrease in systolic blood pressure"
tiff(filename=paste0("b_paternal_w_table_2_r&r.tiff"),
     units="cm", width=45, height=45, res=300)
plot_b=grid.arrange(table_plot_2, table_plot_1, table_plot_3, ncol = 1, bottom=x_title)
plot_b
dev.off()

tiff(filename=paste0("c_paternal_w_table_2_r&r.tiff"),
     units="cm", width=45, height=45, res=300)
plot_c=grid.arrange(table_plot_5, table_plot_8, table_plot_4, ncol = 1, bottom=x_title)
plot_c
dev.off()

# Journal request .eps file

ggsave("plot_a_paternal.eps", plot_a, device="eps",
       units="cm", width=45, height=45)

ggsave("plot_b_paternal.eps", plot_b, device="eps",
       units="cm", width=45, height=45)

ggsave("plot_c_paternal.eps", plot_c, device="eps",
       units="cm", width=45, height=45)


# Results files
results_table_Paternal_DC=as.data.frame(matrix(0,ncol=10, nrow=(nrow(results_all_DC))))
colnames(results_table_Paternal_DC)=c("Outcome","Beta coefficient (95% CI)", "Standard error", "Pvalue", "No. SNPs", "Drug subclass", "Condition", "MR method", "Gene", "rsID")
results_table_Paternal_DC[,"Outcome"]=results_all_DC$Outcome
results_table_Paternal_DC[,"Beta coefficient (95% CI)"]=c(paste0(round(results_all_DC$`Beta coefficient`, 2), " (", round(results_all_DC$`Lower CI`,2), ", ", round(results_all_DC$`Upper CI`,2), ")"))
results_table_Paternal_DC[,"Standard error"]=round(results_all_DC$se,2)
results_table_Paternal_DC[,"MR method"]=results_all_DC$`MR method`
results_table_Paternal_DC[,"Pvalue"]=round(results_all_DC$pval,2)
results_table_Paternal_DC[,"Drug subclass"]=results_all_DC$`Drug subclass`
results_table_Paternal_DC[,"Drug subclass"]=gsub("\n"," ",results_all_DC$`Drug subclass`)
results_table_Paternal_DC[,"No. SNPs"]=results_all_DC$nsnp
results_table_Paternal_DC[,"Gene"]=results_all_DC$Gene
results_table_Paternal_DC[,"rsID"]=results_all_DC$snps
write.csv(results_table_Paternal_DC,"hypertension_subclass_paternal_r&r.csv" )

write.csv(results_all_DC, "gest_r_full_file_pat.csv")

