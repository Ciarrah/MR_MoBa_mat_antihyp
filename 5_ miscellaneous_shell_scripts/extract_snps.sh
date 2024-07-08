dir=""
for i in $dir/*; do eval "j=${i##*/}"; cd $dir/$j; eval "k=*.txt";
for h in $k; do plink2 --bfile .../...(bed/bim/fam location) --extract $j/$h --export A --out output_pathway/output_$j"_"$h ; done; done
