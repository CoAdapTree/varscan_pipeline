### purpose
# use picard to mark/remove duplicates, build bam index for GATK
###

### usage
# 03_mark_build.py /path/to/sortfile /path/to/pooldir/
###

import sys 
import os
from os import path as op
from coadaptree import makedir

thisfile, sortfile, pooldir = sys.argv


# MarkDuplicates
dupdir   = op.join(pooldir,'03_dedup_rg_filtered_indexed_sorted_bamfiles')
pool     = op.basename(pooldir)
samp     = op.basename(sortfile).split("---")[1].split('_R1R2')[0].split(".")[1]
dupfile  = op.join(dupdir,"%s_rd.bam" % samp)
dupflag  = dupfile.replace(".bam",".bam.flagstats")
dupstat  = op.join(dupdir,"%s_rd_dupstat.txt" % samp)

# create sh file
text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=30000M
#SBATCH --ntasks=1
#SBATCH --job-name=mark_%(samp)s
#SBATCH --output=mark_%(samp)s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

# remove dups
module load picard/2.18.9
java -Djava.io.tmpdir=$SLURM_TMPDIR -jar $EBROOTPICARD/picard.jar MarkDuplicates \
I=%(sortfile)s O=%(dupfile)s MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
M=%(dupstat)s REMOVE_DUPLICATES=true
module unload picard

# get more dup stats
module load samtools/1.9
samtools flagstat %(dupfile)s > %(dupflag)s
module unload samtools

# Build bam index for GATK
module load picard/2.18.9
java -jar $EBROOTPICARD/picard.jar BuildBamIndex I=%(dupfile)s
module unload picard

# call next step
source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"

python $HOME/pipeline/04_realignTargetCreator.py %(pooldir)s %(samp)s %(dupfile)s

''' % locals()

# create shdir and file
shdir = op.join(pooldir,'shfiles/03_mark_build_shfiles')
for d in [shdir,dupdir]:
    makedir(d)
file = op.join(shdir,'mark_%(samp)s.sh' % locals())
with open(file,'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(shdir)
print('shdir = ',shdir)
os.system('sbatch %s' % file)

# balance queue
os.system('python $HOME/pipeline/balance_queue.py mark')
os.system('python $HOME/pipeline/balance_queue.py bwa')