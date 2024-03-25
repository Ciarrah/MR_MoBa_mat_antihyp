library(tidyverse);library(stringr);library(readr);library(haven)

genetic_data=read_csv("trio_subset.csv")
phenotypic=read.csv("MBRN_sixmonth_cleaned.csv")

# Load in principal components, merge to genetic file.

pca=read_table2("", col_names = T)
pca$ID_2306[which(pca$Role=="Father")]=gsub("F","",pca$ID_2306[which(pca$Role=="Father")])
pca$ID_2306[which(pca$Role=="Mother")]=gsub("M","",pca$ID_2306[which(pca$Role=="Mother")])

# Split the child id_2306 into barn nr and sentrixid. Split id_2306 into BARN.NR (birth registry file) and preg ID and merge.

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
colnames(Linkage_merge_file)=c("PREG_ID_2306","PAT","MAT","BARN_ID","IID")
pheno_w_link_merge=merge(phenotypic, Linkage_merge_file,by="PREG_ID_2306")

# Tidy environment

rm(FID,MID,child_details)
rm(phenotypic, Linkage_merge_file)
pheno_w_link_merge[,c("MAT","PAT")]=NULL
genetic_data[,c("MAT","PAT")]=NULL

# Merge

geno_pheno_w_link_merge=merge(pheno_w_link_merge, genetic_data, by="IID", all=T)
rm(pheno_w_link_merge, genetic_data)
check=apply(geno_pheno_w_link_merge[,grep("_offspring|_mother|_father",colnames(geno_pheno_w_link_merge))], 1, function(x)all(is.na(x)))
check_2=apply(geno_pheno_w_link_merge[,grep("_offspring|_mother|_father",colnames(geno_pheno_w_link_merge))], 1, function(x)any(is.na(x)))
geno_pheno_w_link_merge=geno_pheno_w_link_merge[which(check==F),]

names(pca)[names(pca) == 'ID_2306']='PREG_ID_2306'
pca_geno_pheno_w_link_merge=merge(geno_pheno_w_link_merge,pca,by="IID")
pca_geno_pheno_w_link_merge=distinct(pca_geno_pheno_w_link_merge)
              
write.table(pca_geno_pheno_w_link_merge,"pheno_geno_pca_merged.txt")
