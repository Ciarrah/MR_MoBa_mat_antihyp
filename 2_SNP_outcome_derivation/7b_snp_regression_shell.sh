#!/bin/bash
#
#SBATCH --job-name=hyp
#SBATCH --output=hyp_output_hyp_%x_%j.out.txt
#SBATCH --error=hyp_error_hyp_%x_%j.out.txt
#SBATCH --account=XXX
#SBATCH --ntasks=XXX
#SBATCH --time=XXX
#SBATCH --mem-per-cpu=XXX

date
cd $PBS_O_WORKDIR
module load R/3.6.0-foss-2019a
Rscript 7a_snp_regression.R $Line
done
