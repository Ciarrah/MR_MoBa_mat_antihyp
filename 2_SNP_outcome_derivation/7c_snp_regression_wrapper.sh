date
File="formatted_hyp_snps.txt"
Lines=$(cat $File)
for Line in $Lines
do
export SNP=$Line
sbatch --export=ALL,Line=$Line 7b_snp_regression_shell.sh
sleep 5
done
