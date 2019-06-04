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
listfile = op.join(aligndir, '%s_realingment_targets.list' % samp)
realbam = op.join(aligndir, '%s_realigned_reads.bam' % samp)

email_text = get_email_info(parentdir, '05')
text = '''#!/bin/bash
#SBATCH --time=7-00:00:00
#SBATCH --mem=8000M
#SBATCH --ntasks=1
#SBATCH --job-name=%(pool)s-%(samp)s-indelRealign
#SBATCH --output=%(pool)s-%(samp)s-indelRealign_%%j.out 
%(email_text)s

module load gatk/3.8
module load java
java -Djava.io.tmpdir=$SLURM_TMPDIR -Xmx8g -jar $EBROOTGATK/GenomeAnalysisTK.jar \
-T IndelRealigner -R %(ref)s -I %(dupfile)s -targetIntervals %(listfile)s -o %(realbam)s
module unload gatk

# sbatch CRISP job if all pooled bamfiles have been created
source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"
export SQUEUE_FORMAT="%%.8i %%.8u %%.12a %%.68j %%.3t %%16S %%.10L %%.5D %%.4C %%.6b %%.7m %%N (%%r)"
python $HOME/pipeline/start_crispANDvarscan.py %(parentdir)s %(pool)s
python $HOME/pipeline/balance_queue.py bedfile

''' % locals()

# create shdir and shfile
shdir = op.join(pooldir, 'shfiles/05_indelRealign_shfiles')
makedir(shdir)
file = op.join(shdir, '%(pool)s-%(samp)s-indelRealign.sh' % locals())
with open(file, 'w') as o:
    o.write("%s" % text)

os.chdir(shdir)
print('shdir = ', shdir)
subprocess.call([shutil.which('sbatch'), file])


balance_queue.main('balance_queue.py', 'indelRealign')
balance_queue.main('balance_queue.py', 'realign')
