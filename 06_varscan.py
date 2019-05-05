"""
### purpose
# use varscan to call snps on individual pool,
# use the GATK to convert vcf to tablefile.txt,
# filter tablefile.txt for non-multiallelic
###

### usage
# python 06_varscan.py /path/to/pooldir/ sampID /p/t/ref.fa /p/t/realigned.bam
###
"""


import sys, os, balance_queue, subprocess, shutil
from os import path as op
from coadaptree import makedir, get_email_info, pklload, fs
from start_crispANDvarscan import sbatch, get_bedfiles


def make_sh(bedfile, shdir, varscandir):
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    bednum = bedfile.split("_")[-1].replace(".bed", "")
    varscanfile = op.join(varscandir, f'{pool}-{samp}-varscan_bedfile_{bednum}.vcf')
    varscanout = varscanfile.replace(".vcf", "_table.txt")
    text = f'''#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --mem=9000M
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --job-name={pool}-{samp}-varscan_bedfile_{bednum}
#SBATCH --output={pool}-{samp}-varscan_bedfile_{bednum}_%j.out 
{email_text}

source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"

# call varscan, convert VCF to table, filter for multi-allelic
module load samtools/1.9
module load java
samtools view -b -L {bedfile} {realbam} > $SLURM_TMPDIR/realigned_{bednum}.bam
samtools mpileup -B -f {ref} $SLURM_TMPDIR/realigned_{bednum}.bam | java -Xmx8g -jar \
$VARSCAN_DIR/VarScan.v2.3.9.jar pileup2cns --min-coverage 8 --p-value 0.05 \
--strand-filter 1 --min-freq-for-hom 0.80 --min-avg-qual 20 --output-vcf 1 > {varscanfile}
module unload samtools

module load gatk/4.1.0.0
gatk VariantsToTable --variant {varscanfile} -F CHROM -F POS -F REF -F ALT -F AF -F QUAL \
-F ADP -F WT -F HET -F HOM -F NC -GF GT -GF GQ -GF SDP -GF DP -GF FREQ -GF PVAL \
-F TYPE -F FILTER -O {varscanout} --split-multi-allelic
module unload gatk

python $HOME/pipeline/balance_queue.py varscan

'''
    file = op.join(shdir, f'{pool}-{samp}-varscan_bedfile_{bednum}.sh')
    with open(file, 'w') as o:
        o.write("%s" % text)
    return file


def create_sh(shdir, varscandir):
    bedfiles = get_bedfiles(parentdir, pool)
    pids = []
    for bedfile in bedfiles:
        file = make_sh(bedfile, shdir, varscandir)
        pids.append(sbatch(file))
    return pids


def create_combine(pids, shdir):
    dependencies = '#SBATCH --dependency=afterok:' + ','.join(pids)
    text = f'''#!/bin/bash
#SBATCH --job-name={pool}-{samp}-combine-varscan
#SBATCH --time=02:59:00
#SBATCH --mem=16000M
#SBATCH --cpus-per-task=1
#SBATCH --output={pool}-{samp}-combine-varscan_%j.out
{dependencies}
{email_text}

source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"
export SQUEUE_FORMAT="%.8i %.8u %.12a %.68j %.3t %16S %.10L %.5D %.4C %.6b %.7m %N (%r)"

python $HOME/pipeline/combine_crispORvarscan.py {pooldir} varscan {samp}

python $HOME/pipeline/balance_queue.py combine-varscan

'''
    file = op.join(shdir, f'{pool}-{samp}-combine-varscan.sh')
    with open(file, 'w') as o:
        o.write("%s" % text)
    sbatch(file)
    print('sbatched combinefile with dependencies: ' + ','.join(pids))


def main():
    # make dirs
    varscandir = op.join(pooldir, 'varscan')
    shdir = op.join(pooldir, 'shfiles/06_varscan_shfiles')
    for d in [varscandir, shdir]:
        makedir(d)

    # get a list of bedfiles, create shfiles, sbatch, return pids
    pids = create_sh(shdir, varscandir)

    # creat file to combine batched varscan calls
    create_combine(pids, shdir)

    # balance queue
    balance_queue.main('balance_queue.py', 'varscan')
    balance_queue.main('balance_queue.py', 'indelRealign')


if __name__ == '__main__':
    thisfile, pooldir, samp, ref, realbam = sys.argv

    email_text = get_email_info(op.dirname(pooldir), '06')
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)

    main()
