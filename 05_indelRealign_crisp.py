"""Create and sbatch gatk indelRealign command files.
Start varscan and crisp if all realigned bamfiles have been made.

### purpose
# use the GATK to realign around indels
###

### usage
# python 05_indelRealign_crisp.py /path/to/pooldir/ sampID
###
"""

import os, sys, balance_queue, subprocess, shutil
from os import path as op
from coadaptree import makedir, get_email_info


thisfile, pooldir, samp, dupfile, ref = sys.argv


# IndelRealigner
pool = op.basename(pooldir)
parentdir = op.dirname(pooldir)
aligndir = op.join(pooldir, '04_realign')
listfile = op.join(aligndir, f'{samp}_realingment_targets.list')
realbam = op.join(aligndir, f'{samp}_realigned_reads.bam')
bash_variables = op.join(parentdir, 'bash_variables')

email_text = get_email_info(parentdir, '05')
text = f'''#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=8000M
#SBATCH --ntasks=1
#SBATCH --job-name={pool}-{samp}-indelRealign
#SBATCH --output={pool}-{samp}-indelRealign_%j.out 
{email_text}

module load gatk/3.8
module load java
export _JAVA_OPTIONS="-Xms256m -Xmx7g"
java -Djava.io.tmpdir=$SLURM_TMPDIR -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T IndelRealigner -R {ref} -I {dupfile} -targetIntervals {listfile} -o {realbam}
module unload gatk

# sbatch CRISP job if all pooled bamfiles have been created
source {bash_variables}
python $HOME/pipeline/start_crispANDvarscan.py {parentdir} {pool}
python $HOME/pipeline/balance_queue.py bedfile {parentdir}

'''

# create shdir and shfile
shdir = op.join(pooldir, 'shfiles/05_indelRealign_shfiles')
makedir(shdir)
file = op.join(shdir, f'{pool}-{samp}-indelRealign.sh')
with open(file, 'w') as o:
    o.write("%s" % text)

os.chdir(shdir)
print('shdir = ', shdir)
subprocess.call([shutil.which('sbatch'), file])


balance_queue.main('balance_queue.py', 'indelRealign', parentdir)
balance_queue.main('balance_queue.py', 'realign', parentdir)
