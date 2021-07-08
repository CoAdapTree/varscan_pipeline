#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=30000M
#SBATCH --ntasks=1
#SBATCH --job-name=DF_pooled-DF_p15-mark
#SBATCH --output=DF_pooled-DF_p15-mark_%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

# remove dups
module load java
module load picard/2.18.9
export _JAVA_OPTIONS="-Xms256m -Xmx27g"
java -Djava.io.tmpdir=$SLURM_TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates I=/scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.004.D709---D505.DF_p15_cap75_kit5_R1R2_trimmed_sorted.bam O=/scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p15_rd.bam MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 M=/scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p15_rd_dupstat.txt REMOVE_DUPLICATES=true

# Build bam index for GATK
java -jar $EBROOTPICARD/picard.jar BuildBamIndex I=/scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p15_rd.bam
module unload picard

# get more dup stats
module load samtools/1.9
samtools flagstat /scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p15_rd.bam > /scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p15_rd.bam.flagstats
module unload samtools

# call next step
source /scratch/lindb/DF_pooled/bash_variables

python $HOME/pipeline/04_realignTargetCreator.py /scratch/lindb/DF_pooled/DF_pooled DF_p15 /scratch/lindb/DF_pooled/DF_pooled/03_dedup_rg_filtered_indexed_sorted_bamfiles/DF_p15_rd.bam

