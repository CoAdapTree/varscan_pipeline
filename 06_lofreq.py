### purpose
# use lofreq to call snps on individual pool,
# use the GATK to convert vcf to tablefile.txt,
# filter tablefile.txt for SNPs
###

### usage
# python 06_lofreq.py /path/to/pooldir/ sampID
###


import sys
import os
from os import path as op
from coadaptree import makedir

thisfile, pooldir, samp, ref, realbam = sys.argv


#LoFreq
lofdir   = op.join(pooldir,'lofreq')
lofile   = op.join(lofdir,'%s_lofreq.vcf.gz' % samp)
lofout   = lofile.replace(".vcf.gz","_table.txt")


text = '''#!/bin/bash
#SBATCH --time=11:59:00
#SBATCH --mem=6000M
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --job-name=lofreq_%(samp)s
#SBATCH --output=lofreq_%(samp)s_%%j.out 
#SBATCH --mail-user=lindb@vcu.edu
#SBATCH --mail-type=FAIL

source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"

# call lofreq, convert VCF to table, filter for multi-allelic
lofreq call-parallel --pp-threads 32 -f %(ref)s -o %(lofile)s %(realbam)s

module load gatk/4.1.0.0
gatk VariantsToTable --variant %(lofile)s -F CHROM -F POS -F REF -F ALT -F AF -F QUAL \
-F DP -F SB -F DP4 -F CONSVAR -F HRUN -F TYPE -F FILTER -O %(lofout)s --split-multi-allelic
module unload gatk

source $HOME/.bashrc
export PYTHONPATH="${PYTHONPATH}:$HOME/pipeline"
python $HOME/pipeline/filter_VariantsToTable.py %(lofout)s

# next step
#python $HOME/pipeline/07_get_snps.py
''' % locals()


shdir = op.join(pooldir,'shfiles/06_lofreq_shfiles')
for d in [lofdir,shdir]:
    makedir(d)
file = op.join(shdir,'lofreq_%(samp)s.sh' % locals())
with open(file,'w') as o:
    o.write("%s" % text)

os.chdir(shdir)
print('shdir = ',shdir)
os.system('sbatch %s' % file)

os.system('python $HOME/pipeline/balance_queue.py lofreq')
os.system('python $HOME/pipeline/balance_queue.py indelRealign')