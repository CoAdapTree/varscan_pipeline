"""
### purpose
# use the GATK to create target intervals for realignment around indels
###

### usage
# python 04_realignTargetCreator.py /path/to/pooldir/ sampID
###
"""

import os, sys, balance_queue, subprocess, shutil
from os import path as op
from coadaptree import makedir, pklload

thisfile, pooldir, samp, dupfile = sys.argv

# RealignerTargetCreator
aligndir = op.join(pooldir, '04_realign')
listfile = op.join(aligndir, '%s_realingment_targets.list' % samp)

# get ref
parentdir = op.dirname(pooldir)
pool = op.basename(pooldir)
ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]

text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=6000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=realign_%(samp)s
#SBATCH --output=realign_%(samp)s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

# realign using the GATK
module load gatk/3.8
java -Djava.io.tmpdir=$SLURM_TMPDIR -Xmx8g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator -R %(ref)s --num_threads 32 -I %(dupfile)s -o %(listfile)s -drf BadMate
module unload gatk

# next step
source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"
python $HOME/pipeline/05_indelRealign_crisp.py %(pooldir)s %(samp)s %(dupfile)s %(ref)s

''' % locals()

# create shdir and shfile
shdir = op.join(pooldir, 'shfiles/04_realignTarget_shfiles')
for d in [aligndir, shdir]:
    makedir(d)
file = op.join(shdir, 'realign_%(samp)s.sh' % locals())
with open(file, 'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(shdir)
print('shdir =', shdir)
subprocess.call([shutil.which('sbatch'), file])

# balance queue
balance_queue.main('balance_queue.py', 'realign')
balance_queue.main('balance_queue.py', 'mark')
