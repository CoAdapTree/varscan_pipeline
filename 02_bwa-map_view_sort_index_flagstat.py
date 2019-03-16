### purpose
# map with bwa, view/sort/index with samtools
###

### usage
# 02_bwa-map_rginfo_mark_realign_lofreq_crisp.py /path/to/ref.fa /path/to/trimmedR1.fastq /path/to/trimmedR2.fastq /sbatch/dir/ sampID
###

### assumes
# outfiles from "bwa index ref.fasta"
#
# export path to lofreq in $HOME/.bashrc
###

import sys, os
from os import path as op
import pickle
from coadaptree import makedir

# get argument inputs
thisfile,ref,r1out,r2out,shdir,samp = sys.argv

### create dirs and filenames
bwashdir  = op.join(shdir,'02_bwa_shfiles')
pooldir   = op.dirname(shdir)
parentdir = op.dirname(pooldir)

# get rginfo
rginfo  = pickle.load(open(op.join(parentdir,'rginfo.pkl'),'rb'))
print('pooldir=',pooldir)
print("RG=",rginfo[samp])
rglb = rginfo[samp]['rglb']
rgpl = rginfo[samp]['rgpl']
rgsm = rginfo[samp]['rgsm']

#bwa: fastq -> sam
sam      = op.basename(r1out).replace("R1_trimmed.fastq.gz","R1R2_trimmed.sam")
samdir   = op.join(pooldir,'02a_samfiles')
samfile  = op.join(samdir,sam)
#samtools view: sam -> bam
bam      = op.basename(samfile).replace('.sam','.bam')
bamdir   = op.join(pooldir,'02b_bamfiles')
bamfile  = op.join(bamdir,bam)
#samtools sort: bamfile -> sortfile
sort     = op.basename(bamfile).replace('.bam','_sorted.bam')
sortdir  = op.join(pooldir,'02c_sorted_bamfiles')
sortfile = op.join(sortdir,sort)
flagfile = op.join(sortdir,sort.replace('.bam','.bam.flagstats'))


for d in [bwashdir,samdir,bamdir,sortdir]:
    makedir(d)

# send it off
text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=30000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=bwa_%(samp)s
#SBATCH --output=bwa_%(samp)s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

# get RGID and RGPU
RGID=$(zcat %(r1out)s | head -n1 | sed 's/:/_/g' | cut -d "_" -f1,2,3,4)
RGPU=$RGID.%(rglb)s

# map, sam to bam, sort by coordinate, index
module load bwa/0.7.17
bwa mem -t 32 -M -R "@RG\\tID:$RGID\\tSM:%(rgsm)s\\tPL:%(rgpl)s\\tLB:%(rglb)s\\tPU:$RGPU" \
%(ref)s %(r1out)s %(r2out)s > %(samfile)s
module unload bwa

module load samtools/1.9
samtools view -@ 32 -q 20 -F 0x0004 -Sb %(samfile)s > %(bamfile)s
samtools sort -@ 32 %(bamfile)s > %(sortfile)s
samtools index %(sortfile)s
samtools flagstat %(sortfile)s > %(flagfile)s
module unload samtools

# mark and build
source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"
python $HOME/pipeline/03_mark_build.py %(sortfile)s %(pooldir)s 
''' % locals()

# create shfile
qsubfile = op.join(bwashdir,'bwa_%s.sh' % samp) 
with open(qsubfile,'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(bwashdir)
print('shdir = ',shdir)
os.system("sbatch %s" % qsubfile)

os.system('python $HOME/pipeline/balance_queue.py bwa')
os.system('python $HOME/pipeline/balance_queue.py trim')