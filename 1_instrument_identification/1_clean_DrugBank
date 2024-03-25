library(readr);library(tidyr)
setwd()

pharmacological_file_pathway="pharmacologically_active.csv"
vocab_file_pathway="drugbank_vocabulary.csv"

format_drugbank=function(vocab_file_pathway, pharmacological_file_pathway)
{
  vocab=read_csv(vocab_file_pathway)
  vocab$substance=paste0(vocab$`Common name`,"|",vocab$Synonyms)
  vocab=vocab[,c("DrugBank ID","substance")]
  vocab$substance=gsub("\\|",";",vocab$substance) 
  vocab=separate_rows(vocab,substance,sep = ";")
  
  pharma=read_csv(pharmacological_file_pathway)
  pharma=pharma[,c("Gene Name","Drug IDs")]
  colnames(pharma)=c("gene","DrugBank ID")
  pharma=separate_rows(pharma, `DrugBank ID`)
  
  complete_file=merge(vocab,pharma,by=c("DrugBank ID"),all.x = TRUE)
  complete_file=complete_file[complete_file$substance!="NA",]
  complete_file$substance=tolower(complete_file$substance)
  complete_file=unique(complete_file)
  
  return(complete_file)
}
