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
from coadaptree import makedir, pklload, get_email_info

thisfile, pooldir, samp, dupfile = sys.argv

# RealignerTargetCreator
aligndir = op.join(pooldir, '04_realign')
listfile = op.join(aligndir, '%s_realingment_targets.list' % samp)

# get ref
parentdir = op.dirname(pooldir)
pool = op.basename(pooldir)
ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]

email_text = get_email_info(parentdir, '04')
text = '''#!/bin/bash
#SBATCH --time=3-00:00:00
#SBATCH --mem=8000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=%(pool)s-%(samp)s-realign
#SBATCH --output=%(pool)s-%(samp)s-realign_%%j.out 
%(email_text)s

# realign using the GATK
module load gatk/3.8
java -Djava.io.tmpdir=$SLURM_TMPDIR -Xmx8g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator -R %(ref)s --num_threads 32 -I %(dupfile)s -o %(listfile)s
module unload gatk

# next step
source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"
export SQUEUE_FORMAT="%%.8i %%.8u %%.12a %%.68j %%.3t %%16S %%.10L %%.5D %%.4C %%.6b %%.7m %%N (%%r)"
python $HOME/pipeline/05_indelRealign_crisp.py %(pooldir)s %(samp)s %(dupfile)s %(ref)s

''' % locals()

# create shdir and shfile
shdir = op.join(pooldir, 'shfiles/04_realignTarget_shfiles')
for d in [aligndir, shdir]:
    makedir(d)
file = op.join(shdir, '%(pool)s-%(samp)s-realign.sh' % locals())
with open(file, 'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(shdir)
print('shdir =', shdir)
subprocess.call([shutil.which('sbatch'), file])

# balance queue
balance_queue.main('balance_queue.py', 'realign')
balance_queue.main('balance_queue.py', 'mark')
