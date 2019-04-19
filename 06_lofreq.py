"""
### purpose
# use lofreq to call snps on individual pool,
# use the GATK to convert vcf to tablefile.txt,
# filter tablefile.txt for SNPs
###

### usage
# python 06_lofreq.py /path/to/pooldir/ sampID
###
"""


import sys, os, balance_queue, subprocess, shutil
from os import path as op
from coadaptree import makedir, get_email_info, pklload, fs
from start_crisp import sbatch, get_bedfiles


def make_sh(bedfile, shdir, lofdir):
    bednum = bedfile.split("_")[-1].replace(".bed", "")
    lofile = op.join(lofdir, f'{pool}-{samp}-lofreq_bedfile_{bednum}.vcf.gz')
    lofout = lofile.replace(".vcf.gz", "_table.txt")
    text = f'''#!/bin/bash
#SBATCH --time=23:59:00
#SBATCH --mem=9000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name={pool}-{samp}-lofreq_bedfile_{bednum}
#SBATCH --output={pool}-{samp}-lofreq_bedfile_{bednum}_%j.out 
{email_text}

source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"

# call lofreq, convert VCF to table, filter for multi-allelic
lofreq call-parallel --pp-threads 32 -f {ref} -o {lofile} {realbam} -l {bedfile}

module load gatk/4.1.0.0
gatk VariantsToTable --variant {lofile} -F CHROM -F POS -F REF -F ALT -F AF -F QUAL \
-F DP -F SB -F DP4 -F CONSVAR -F HRUN -F TYPE -F FILTER -O {lofout} --split-multi-allelic
module unload gatk

python $HOME/pipeline/balance_queue.py lofreq

'''
    file = op.join(shdir, f'{pool}-{samp}-lofreq_bedfile_{bednum}.sh')
    with open(file, 'w') as o:
        o.write("%s" % text)
    return file


def create_sh(shdir, lofdir):
    bedfiles = get_bedfiles(parentdir, pool)
    pids = []
    for bedfile in bedfiles:
        file = make_sh(bedfile, shdir, lofdir)
        pids.append(sbatch(file))
    return pids


def create_combine(pids, shdir):
    dependencies = '#SBATCH --dependency=afterok:' + ','.join(pids)
    text = f'''#!/bin/bash
#SBATCH --job-name={pool}-{samp}-combine-lofreq
#SBATCH --time=02:59:00
#SBATCH --mem=16000M
#SBATCH --cpus-per-task=1
#SBATCH --output={pool}-{samp}-combine-lofreq_%j.out
{dependencies}
{email_text}

source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"
export SQUEUE_FORMAT="%.8i %.8u %.12a %.68j %.3t %16S %.10L %.5D %.4C %.6b %.7m %N (%r)"

python $HOME/pipeline/combine_crispORlofreq.py {pooldir} lofreq {samp}

python $HOME/pipeline/balance_queue.py combine-lofreq

'''
    file = op.join(shdir, f'{pool}-{samp}-combine-lofreq.sh')
    with open(file, 'w') as o:
        o.write("%s" % text)
    sbatch(file)
    print('sbatched combinefile with dependencies: ' + ','.join(pids))


def main():
    # make dirs
    lofdir = op.join(pooldir, 'lofreq')
    shdir = op.join(pooldir, 'shfiles/06_lofreq_shfiles')
    for d in [lofdir, shdir]:
        makedir(d)

    # get a list of bedfiles, create shfiles, sbatch, return pids
    pids = create_sh(shdir, lofdir)

    # creat file to combine batched lofreq calls
    create_combine(pids, shdir)

    # balance queue
    balance_queue.main('balance_queue.py', 'lofreq')
    balance_queue.main('balance_queue.py', 'indelRealign')


if __name__ == '__main__':
    thisfile, pooldir, samp, ref, realbam = sys.argv

    email_text = get_email_info(op.dirname(pooldir), '06')
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)

    main()
