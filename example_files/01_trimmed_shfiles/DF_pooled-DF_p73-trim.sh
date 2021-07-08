#!/bin/bash
#SBATCH --job-name=DF_pooled-DF_p73-trim
#SBATCH --time=02:59:00
#SBATCH --mem=5000M
#SBATCH --cpus-per-task=16
#SBATCH --output=DF_pooled-DF_p73-trim_%j.out
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source /scratch/lindb/DF_pooled/bash_variables

module load fastp/0.19.5

fastp -i /scratch/lindb/DF_pooled/DF_pooled/NS.1195.003.D709---D507.DF_p73_cap57_kit5_R1.fastq.gz -o /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.003.D709---D507.DF_p73_cap57_kit5_R1_trimmed.fastq.gz -I /scratch/lindb/DF_pooled/DF_pooled/NS.1195.003.D709---D507.DF_p73_cap57_kit5_R2.fastq.gz -O /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.003.D709---D507.DF_p73_cap57_kit5_R2_trimmed.fastq.gz -g --cut_window_size 5 --cut_mean_quality 30 --n_base_limit 20 --length_required 75 -h /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.003.D709---D507.DF_p73_cap57_kit5__trimmed_R1_R2_stats.html --cut_by_quality3 --thread 16 --json /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.003.D709---D507.DF_p73_cap57_kit5__trimmed_R1_R2.json --adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC --adapter_sequence_r2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT > /scratch/lindb/DF_pooled/DF_pooled/01_trimmed/NS.1195.003.D709---D507.DF_p73_cap57_kit5__trimmed_R1_R2_stats.log

# once finished, map using bwa mem
python $HOME/pipeline/02_bwa-map_view_sort_index_flagstat.py /scratch/lindb/DF_pooled DF_p73

