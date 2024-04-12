#!/bin/bash
#
#SBATCH --job-name=r_r_hyp
#SBATCH --output=output_r_r_%x_%j.out.txt
#SBATCH --error=error_hyp_r_r_%x_%j.out.txt
#SBATCH --account=XXX
#SBATCH --ntasks=1
#SBATCH --time=03:00:00
#SBATCH --mem-per-cpu=250G

date
cd $PBS_O_WORKDIR
module load R/3.6.0-foss-2019a
Rscript link_pheno_geno_hyp.R
done
