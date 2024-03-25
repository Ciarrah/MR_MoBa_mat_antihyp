library(data.table);library(ggplot2);library(stringr);library(tidyverse);library(haven)
`%!in%` = Negate(`%in%`)

output_snps=fread("")               # .raw output file from PLINK extract by gene chr_pos
consent=read.delim("")

output_snps=output_snps[which(output_snps$IID%in%consent$IID),]

output_snps=output_snps[,-which(colnames(output_snps)=="PHENOTYPE")]
snp_names=grep("rs", colnames(output_snps), value=T)
other_cols=colnames(output_snps)[-which(colnames(output_snps)%in%c(snp_names,"FID"))]

fathers=output_snps[which(output_snps$SEX==1&output_snps$PAT==0),]
mothers=output_snps[which(output_snps$SEX==2&output_snps$MAT==0),]
offspring=output_snps[which(output_snps$PAT!=0|output_snps$MAT!=0),]

colnames(fathers)[which(colnames(fathers)%in%snp_names)]=paste0(colnames(fathers)[which(colnames(fathers)%in%snp_names)],"_fathers")
colnames(mothers)[which(colnames(mothers)%in%snp_names)]=paste0(colnames(mothers)[which(colnames(mothers)%in%snp_names)],"_mothers")
colnames(offspring)[which(colnames(offspring)%in%snp_names)]=paste0(colnames(offspring)[which(colnames(offspring)%in%snp_names)],"_offspring")
mothers[,c("PAT","MAT","SEX")]=fathers[,c("PAT","MAT","SEX")]=NULL
new_subset=merge(mothers,fathers,by="FID")

# Father id is y mother id is x
colnames(new_subset)[which(colnames(new_subset)=="IID.y")]="PAT"
colnames(new_subset)[which(colnames(new_subset)=="IID.x")]="MAT"

# Some family IDs have differing fathers

trio_subset=NULL
for (i in 1:nrow(offspring))
{
  j=which(new_subset$MAT==offspring$MAT[i]&new_subset$PAT==offspring$PAT[i])
  single_row=merge(offspring[i,],new_subset[j,])
  trio_subset=rbind(trio_subset, single_row)
  
  # For those with father only
  
  if (offspring$MAT[i]==0&offspring$PAT[i]!=0)
  {
    j=which(new_subset$FID==offspring$FID[i]&new_subset$PAT==offspring$PAT[i])
    single_row=merge(offspring[i,],new_subset[j,])
    trio_subset=rbind(trio_subset, single_row)
  }
  
  # For those with mother only
  
  else if (offspring$PAT[i]==0&offspring$MAT[i]!=0)
  {
    j=which(new_subset$FID==offspring$FID[i]&new_subset$MAT==offspring$MAT[i])
    single_row=merge(offspring[i,],new_subset[j,])
    trio_subset=rbind(trio_subset, single_row)
  }
  print(i) 
}
write.csv(trio_subset,"trio_subset.csv")
