#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=30000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=DF_pooled-DF_p44-realign
#SBATCH --output=DF_pooled-DF_p44-realign_%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

# realign using the GATK
module load java
module load gatk/3.8
export _JAVA_OPTIONS="-Xms256m -Xmx28g"
java -Djava.io.tmpdir=$SLURM_TMPDIR -jar $EBROOTGATK/GenomeAnalysisTK.jar -T RealignerTargetCreator -R /scratch/lindb/DF_2018_12_edit/DF_ref_edit.fasta --num_threads 32 -I /scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p44_rd.bam -o /scratch/lindb/DF_pooled/DF_pooled/04_realign/DF_p44_realingment_targets.list
module unload gatk

# next step
source /scratch/lindb/DF_pooled/bash_variables
python $HOME/pipeline/05_indelRealign.py /scratch/lindb/DF_pooled/DF_pooled DF_p44 /scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p44_rd.bam /scratch/lindb/DF_2018_12_edit/DF_ref_edit.fasta

