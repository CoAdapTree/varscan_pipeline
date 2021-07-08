#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=8000M
#SBATCH --ntasks=1
#SBATCH --job-name=DF_pooled-DF_p43-indelRealign
#SBATCH --output=DF_pooled-DF_p43-indelRealign_%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

module load java
module load gatk/3.8
export _JAVA_OPTIONS="-Xms256m -Xmx7g"
java -Djava.io.tmpdir=$SLURM_TMPDIR -jar $EBROOTGATK/GenomeAnalysisTK.jar -T IndelRealigner -R /scratch/lindb/DF_2018_12_edit/DF_ref_edit.fasta -I /scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p43_rd.bam -targetIntervals /scratch/lindb/DF_pooled/DF_pooled/04_realign/DF_p43_realingment_targets.list -o /scratch/lindb/DF_pooled/DF_pooled/04_realign/DF_p43_realigned_reads.bam
module unload gatk

# sbatch varscan jobs if all pooled bamfiles have been created
source /scratch/lindb/DF_pooled/bash_variables
python $HOME/pipeline/start_varscan.py /scratch/lindb/DF_pooled DF_pooled
python $HOME/pipeline/balance_queue.py bedfile /scratch/lindb/DF_pooled

