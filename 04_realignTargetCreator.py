"""Create and sbatch gatk realignTargetCreator command files.

### purpose
# use the GATK to create target intervals for realignment around indels
###

### usage
# python 04_realignTargetCreator.py /path/to/pooldir/ sampID
###
"""

import os, sys, subprocess, shutil
from os import path as op
from coadaptree import makedir, pklload, get_email_info

thisfile, pooldir, samp, dupfile = sys.argv

# RealignerTargetCreator
aligndir = op.join(pooldir, '04_realign')
listfile = op.join(aligndir, f'{samp}_realingment_targets.list')

# get ref
parentdir = op.dirname(pooldir)
pool = op.basename(pooldir)
ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
bash_variables = op.join(parentdir, 'bash_variables')

email_text = get_email_info(parentdir, '04')
text = f'''#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=30000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name={pool}-{samp}-realign
#SBATCH --output={pool}-{samp}-realign_%j.out 
{email_text}

# realign using the GATK
module load gatk/3.8
module load java/1.8.0_192
export _JAVA_OPTIONS="-Xms256m -Xmx28g"
java -Djava.io.tmpdir=$SLURM_TMPDIR -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T RealignerTargetCreator -R {ref} --num_threads 32 -I {dupfile} -o {listfile}
module unload gatk

# next step
source {bash_variables}
python $HOME/pipeline/05_indelRealign.py {pooldir} {samp} {dupfile} {ref}

'''

# create shdir and shfile
shdir = op.join(pooldir, 'shfiles/04_realignTarget_shfiles')
for d in [aligndir, shdir]:
    makedir(d)
file = op.join(shdir, f'{pool}-{samp}-realign.sh')
with open(file, 'w') as o:
    o.write("%s" % text)

# sbatch file
os.chdir(shdir)
print('shdir =', shdir)
subprocess.call([shutil.which('sbatch'), file])

# balance queue
balance_queue = op.join(os.environ['HOME'], 'pipeline/balance_queue.py')
subprocess.call([sys.executable, balance_queue, 'realign', parentdir])
subprocess.call([sys.executable, balance_queue, 'mark', parentdir])
