library(data.table);library(ggplot2);library(readr)
`%!in%` = Negate(`%in%`)

pheno_geno=fread("pheno_geno_pca_merged.txt")
pheno_geno=as.data.frame(pheno_geno)

# Calculate the mean allele frequency for each SNP and divide by two.

results=NULL
results=colMeans(pheno_geno[,grep("_offspring|_mothers|_fathers",colnames(pheno_geno),value=T)],na.rm = TRUE)
results=results/2
results=as.data.frame(results)
colnames(results)="mean_allele_frq"
results$sample.size=lapply(lapply(pheno_geno[grep("_offspring|_mothers|_fathers",colnames(pheno_geno))],table),sum)
results$SNP=rownames(results)
results$sample.size=as.character(results$sample.size)
results=results[,c("SNP", "mean_allele_frq","sample.size")]
write.csv(results,file="trio_allele_frequencies.csv")

off_freq=as.data.frame(results[grep("offspring",rownames(results)),])
mat_freq=as.data.frame(results[grep("mothers",rownames(results)),])
pat_freq=as.data.frame(results[grep("fathers",rownames(results)),])

rownames(off_freq)=gsub("_.*","",grep("offspring",rownames(results),value=T))
rownames(mat_freq)=gsub("_.*","",grep("mothers",rownames(results),value=T))
rownames(pat_freq)=gsub("_.*","",grep("fathers",rownames(results),value=T))

off_freq$SNP=rownames(off_freq)
mat_freq$SNP=rownames(mat_freq)
pat_freq$SNP=rownames(pat_freq)

colnames(off_freq)=c("SNP","Off_AF","Off_ss")
colnames(mat_freq)=c("SNP","Mat_AF","Mat_ss")
colnames(pat_freq)=c("SNP","Pat_AF","Pat_ss")

DF=merge(merge(pat_freq,off_freq,by="SNP"),mat_freq,by="SNP")
fwrite(DF, "trio_allele_freq.csv")

# Compare the frequencies between pairs using scatter plots.

pdf("mean_trio_allele_comparison.pdf")
ggplot(DF,aes(x=Off_AF,y=Mat_AF))+geom_point()+labs(x="Offspring mean allele frequency", y="Maternal allele frequency")
ggplot(DF,aes(x=Off_AF,y=Pat_AF))+geom_point()+labs(x="Offspring mean allele frequency", y="Paternal allele frequency")
ggplot(DF,aes(x=Pat_AF,y=Mat_AF))+geom_point()+labs(x="Paternal mean allele frequency", y="Maternal allele frequency")
dev.off()

# Allele frequencies appear concurrent.

## Load the 1kg data and compare frequencies, majority of snps available, manually check on SNPdb if not.

output.afreq=fread("")
colnames(output.afreq)=c("#CHROM","SNP",	"REF","ALT","ALT_FREQS","OBS_CT")
output.afreq=output.afreq[which(output.afreq$SNP%in%DF$SNP),]
length(which(output.afreq$SNP%in%DF$SNP))

# Check the direction of effect and reference allele for MoBa and 1kg are the same.

moba_snp_direction=unique(sub("mothers|fathers|offspring","",colnames(pheno_geno)))
moba_snp_direction=as.data.frame(grep("_A_|_G_|_C_|_T_",moba_snp_direction,value=T))
colnames(moba_snp_direction)="SNP"
moba_snp_direction$effect_allele=NA
moba_snp_direction$effect_allele=sapply(strsplit(moba_snp_direction$SNP,"_"),"[",2)
moba_snp_direction$SNP=gsub("_A_|_C_|_G_|_T_","",moba_snp_direction$SNP)

write_csv(moba_snp_direction, "moba_snps_effect_alleles.csv")

moba_snp_direction=moba_snp_direction[which(moba_snp_direction$SNP%in%output.afreq$SNP),]
moba_snp_direction=moba_snp_direction[,c("SNP","effect_allele")]

DF=merge(DF, moba_snp_direction, by="SNP")
DF=merge(DF, output.afreq, by="SNP")
file=DF[which(DF$effect_allele!=DF$REF),]
file_2=DF[which(DF$effect_allele==DF$REF),]
file_2$ALT_FREQS=1-file_2$ALT_FREQS
DF=rbind(file, file_2)

colnames(DF)[which(colnames(DF)=="ALT_FREQS")]="aligned_freq"

# XXX SNPs not aligned the same way

# Create the aligned eaf between 1kg and moba, want x=y trend

pdf("mean_trio_1kg_allele_comparison.pdf", onefile=T)
ggplot(DF,aes(x=Off_AF,y=aligned_freq))+geom_point()+labs(x="Offspring mean allele frequency", y="1000 genomes mean allele frequency")
ggplot(DF,aes(x=Pat_AF,y=aligned_freq))+geom_point()+labs(x="Paternal mean allele frequency", y="1000 genomes mean allele frequency")
ggplot(DF,aes(x=Mat_AF,y=aligned_freq))+geom_point()+labs(x="Maternal mean allele frequency", y="1000 genomes mean allele frequency")
dev.off()

# General linear trend, expected.

# Create summary statistics file for genetic MoBa data

sum_tab=as.data.frame(matrix(NA,nrow=nrow(DF),ncol=8))
colnames(sum_tab)=c("RSID","CHR","POS","Effect allele","Other allele","Offspring effect allele",
                    "Maternal effect allele", "Paternal effect allele")

# Add in MoBa mean eaf data

sum_tab[,c("RSID","Offspring effect allele","Maternal effect allele",
           "Paternal effect allele","Effect allele","Other allele")]=
  DF[,c("SNP","Off_AF","Mat_AF",'Pat_AF', "effect_allele","ALT")]

for (i in file$SNP)
{
  sum_tab[which(sum_tab$RSID==i),"Effect allele"]=file$ALT[which(file$SNP==i)]
  sum_tab[which(sum_tab$RSID==i),"Other allele"]=file$REF[which(file$SNP==i)]
  print(which(file$SNP==i))
}

# Check order is the same

table((order(sum_tab$RSID)==order(DF$SNP)))

for (k in output.afreq$SNP)
{
  sum_tab$CHR[which(sum_tab$RSID==k)]=output.afreq[which(output.afreq$SNP==k),"#CHROM"]
}

bim_file=fread("")
colnames(bim_file)=c("CHR","RSID","POS","BP_COR","REF_AL","ALT_AL")

# Restrict to RSID of interest.

bim_file=bim_file[which(bim_file$RSID%in%sum_tab$RSID),c("RSID","BP_COR")]
sum_tab=merge(sum_tab,bim_file,by="RSID")
sum_tab$POS=sum_tab$BP_COR
sum_tab$BP_COR=NULL
#sum_tab=apply(sum_tab,2,as.character)

## Add back mat, pat, off freq for rsid missing in 1kg

additional_rows_O=off_freq[which(off_freq$SNP%!in%output.afreq$SNP),c("SNP","Off_AF")]
additional_rows_m=mat_freq[which(mat_freq$SNP%!in%output.afreq$SNP),c("SNP","Mat_AF")]
additional_rows_f=pat_freq[which(pat_freq$SNP%!in%output.afreq$SNP),c("SNP","Pat_AF")]
table(order(additional_rows_f$SNP)==order(additional_rows_O$SNP))
table(order(additional_rows_m$SNP)==order(additional_rows_O$SNP))

moba_snp_direction=unique(sub("mothers|fathers|offspring","",colnames(pheno_geno)))
moba_snp_direction=as.data.frame(grep("_A_|_G_|_C_|_T_",moba_snp_direction,value=T))
colnames(moba_snp_direction)="SNP"
moba_snp_direction$effect_allele=NA
moba_snp_direction$effect_allele=sapply(strsplit(moba_snp_direction$SNP,"_"),"[",2)
moba_snp_direction$SNP=gsub("_A_|_C_|_G_|_T_","",moba_snp_direction$SNP)

additional_rows=as.data.frame(cbind(additional_rows_f$SNP,NA,NA,
                      moba_snp_direction[which(moba_snp_direction$SNP%in%additional_rows_f$SNP),"effect_allele"],
                      NA,additional_rows_O$Off_AF, additional_rows_m$Mat_AF, additional_rows_f$Pat_AF))

colnames(additional_rows)=colnames(sum_tab)

sum_tab=rbind(sum_tab, additional_rows)
df=apply(sum_tab,2,as.character)
write.csv(df,"summary_table_genetic_trio.csv")
