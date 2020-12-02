#!/bin/bash
#SBATCH --job-name=DF_pooled-combine-varscan
#SBATCH --time=12:00:00
#SBATCH --mem=200000M
#SBATCH --cpus-per-task=1
#SBATCH --output=DF_pooled-combine-varscan_%j.out
#SBATCH --dependency=afterok:25170544,25170545,25170546,25170547,25170548,25170549,25170550,25170551,25170552
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL


source /scratch/lindb/DF_pooled/bash_variables

python $HOME/pipeline/combine_varscan.py /scratch/lindb/DF_pooled/DF_pooled varscan DF_pooled

