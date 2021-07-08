#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --mem=55000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=DF_pooled-DF_p7-bwa
#SBATCH --output=DF_pooled-DF_p7-bwa_%j.out
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

# get RGID and RGPU
RGID=$(zcat /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1_trimmed.fastq.gz | head -n1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
RGPU=$RGID.DF.cap18.kit3

# map, sam to bam, sort by coordinate, index
module load bwa/0.7.17
bwa mem -t 32 -M -R "@RG\tID:$RGID\tSM:DF.p7.cap18.kit3\tPL:ILLUMINA\tLB:DF.cap18.kit3\tPU:$RGPU" /scratch/lindb/DF_2018_12_edit/DF_ref_edit.fasta /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1_trimmed.fastq.gz /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R2_trimmed.fastq.gz > /scratch/lindb/DF_pooled/DF_pooled/02a_samfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed.sam
module unload bwa

module load samtools/1.9
samtools view -@ 32 -q 20 -F 0x0004 -f 0x0002 -Sb /scratch/lindb/DF_pooled/DF_pooled/02a_samfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed.sam > /scratch/lindb/DF_pooled/DF_pooled/02b_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed.bam
samtools sort -@ 32 /scratch/lindb/DF_pooled/DF_pooled/02b_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed.bam > /scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed_sorted.bam
samtools index /scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed_sorted.bam
samtools flagstat /scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed_sorted.bam > /scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed_sorted.bam.flagstats
module unload samtools

module load bedtools/2.27.1
bedtools bamtobed -i /scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed_sorted.bam > /scratch/lindb/DF_pooled/DF_pooled/02c_sorted_bamfiles/NS.1195.001.D703---D507.DF_p7_cap18_kit3_R1R2_trimmed_sortedbam.coord
module unload bedtools



# mark and build
source /scratch/lindb/DF_pooled/bash_variables
python $HOME/pipeline/03_mark_build.py /scratch/lindb/DF_pooled/DF_pooled DF_p7
