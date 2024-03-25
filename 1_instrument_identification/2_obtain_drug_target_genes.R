library(readxl);library(tidyr);library(dplyr)
`%!in%` = Negate(`%in%`)
setwd()

# Format BNF file. Remove dosage, retain only drug substance. Ensure drug substance is compitable with
# complete file to merge. Remove polypharmacy.

# dm+d antihypertensive search for openprescribing_search_hypertension.csv :
# https://openprescribing.net/dmd/advanced-search/vmp/?search=%5B%22or%22%2C%5B%5B%22bnf_code%22%2C%22begins_with%22%2C%2202040%22%5D%2C%5B%22bnf_
# code%22%2C%22begins_with%22%2C%22020203%22%5D%2C%5B%22bnf_code%22%2C%22begins_with%22%2C%22020602%22%5D%2C%5B%22bnf_code%22%2C%22begins_with%22%2C%22020201%22
# %5D%2C%5B%22bnf_code%22%2C%22begins_with%22%2C%22020202%22%5D%2C%5B%22bnf_code%22%2C%22begins_with%22%2C%2202050%22%5D%2C%5B%22bnf_code%22%2C%22begins_with%22
# %2C%220401020R0%22%5D%2C%5B%22nm%22%2C%22contains%22%2C%22methyclothiazide%22%5D%2C%5B%22nm%22%2C%22contains%22%2C%22hydroflumethiazide%22%5D%2C%5B%22nm%22%2C%22
# contains%22%2C%22metirosine%22%5D%2C%5B%22nm%22%2C%22contains%22%2C%22mibefradul%22%5D%2C%5B%22nm%22%2C%22contains%22%2C%22betanidine%22%5D%2C%5B%22nm%22%2C%22
# contains%22%2C%22perhexiline%22%5D%2C%5B%22bnf_code%22%2C%22begins_with%22%2C%220902011U0%22%5D%2C%5B%22nm%22%2C%22contains%22%2C%22tadalafil%22%5D%2C%5B%22
# nm%22%2C%22contains%22%2C%22tamsulosin%22%5D%5D%5D&include=unavailable&include=invalid&include=no_bnf_code

# Clean and format DrugBank file

source("1_clean_DrugBank.R")
pharmacological_file_pathway="pharmacologically_active.csv"
vocab_file_pathway="drugbank_vocabulary.csv"
complete_file=format_drugbank(vocab_file_pathway, pharmacological_file_pathway)

# Format BNF file

BNF_drug_classes=read.csv("openprescribing_search_hypertension.csv")
polypharmacy=BNF_drug_classes[grep("\\b\\S+\\s/\\s\\S+\\b", BNF_drug_classes$name),]
BNF_drug_classes=BNF_drug_classes[-which(BNF_drug_classes$name%in%polypharmacy$name),]

# BNF drug name formating, remove dosage and make lowercase.
remove_dose=function(x) {
  numeric_index=regexpr("[0-9]", x)
  if (numeric_index>0) {
    return(substr(x, 1, numeric_index-1))
  } else {
    return(x)
  }
}
BNF_drug_classes$name=sapply(BNF_drug_classes$name, remove_dose)
BNF_drug_classes$name=tolower(BNF_drug_classes$name)

# Remove excess characters in BNF chapter (pertaining to formulation etc),remove duplicates and replace capitals to match complete_file format.

BNF_drug_classes$index=regexpr("[A-Za-z]", BNF_drug_classes$bnf_code)
BNF_drug_classes$bnf_code=substr(BNF_drug_classes$bnf_code, 1, BNF_drug_classes$index-1)
BNF_drug_classes[,c("vmp_id","invalid","index")]=NULL
BNF_drug_classes=distinct(BNF_drug_classes)

# 106 unique substances found in d+dm 
colnames(BNF_drug_classes)=tolower(colnames(BNF_drug_classes))

# Remove trailing/leading whitespace in both files

BNF_drug_classes$name=trimws(BNF_drug_classes$name, "l")
BNF_drug_classes$name=trimws(BNF_drug_classes$name, "r")

complete_file$substance=trimws(complete_file$substance, "l")
complete_file$substance=trimws(complete_file$substance, "r")

length(which(BNF_drug_classes$name%in%complete_file$substance))

# 92 direct matches, review non matches first

# Merge file to obtain genes of interest.

test_file=merge(complete_file, BNF_drug_classes, by.x="substance", by.y="name", all.y = T)

# Manually review names for substances that were not linked.

length(which(BNF_drug_classes$name%!in%complete_file$substance))
manually_check=BNF_drug_classes[which(BNF_drug_classes$name%!in%complete_file$substance),]

# Remove "generic sevikar" - already included as Olmesartan medoxomil. Remove "generic titrace" - already included
# as ramipril. Remove "generic kalten" and "generic moducren" - included as amiloride.

manually_check=manually_check[-grep("generic",manually_check$name),]
manually_check[grep("perindopril",manually_check$name),"name"]="perindopril"
manually_check=unique(manually_check)
manually_check[grep("olmesartan",manually_check$name),"name"]="olmesartan"
manually_check[grep("potassium", manually_check$name),"name"]="potassium chloride"

# Remove, not in complete file nor mentioned in VW codelist. Extremely rarely used substance, emergencies only.

manually_check=manually_check[-grep("trimetaphan camsilate",manually_check$name),]

# Co-tenidone is a combination of atenolol and chlortalidone, remove. Co-prenozide is Oxprenolol hydrochloride/cyclopenthiazide), remove
# Co-zidocapt is captopril and hydrochlorothiazide, remove.

manually_check=manually_check[-grep("co-", manually_check$name),]
test_2=merge(complete_file, manually_check, by.x="substance", by.y="name", all.y=T)

# Checks complete, bind two files. Retain only rows with genetic target.

merged_file=rbind(test_file, test_2)
length(unique(merged_file$substance))

# 87 unique drug substances

merged_file=subset(merged_file, complete.cases(gene))
merged_file=unique(merged_file)

# 87 unique substances across 73 genes

length(unique(merged_file$substance))
length(unique(merged_file$gene))

merged_file$condition="Hypertension"

# Restrict to DrugBank indicated genes

write_csv(merged_file,"identified_genes_hypertension_190523.csv")

