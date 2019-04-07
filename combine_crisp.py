"""
### purpose
# combine the bedfiles from parallelized crisp runs
###

### usage
# python combine_crisp.py /path/to/pooldir/
###
"""

import os, sys, time, random, pandas as pd
from os import path as op
from coadaptree import fs, pklload
from filter_VariantsToTable import main as remove_multiallelic
from start_crisp import getfiles


def checkjobs(pooldir):
    pool = op.basename(pooldir)
    samps = pklload(op.join(op.dirname(pooldir), 'poolsamps.pkl'))[pool]
    shdir = op.join(pooldir, 'shfiles/crisp')
    files = getfiles(samps, shdir, 'crisp_bedfile')
#     check_queue(files.values(), pooldir)  # make sure job isn't in the queue (running or pending)
#     check_seff(files.values())  # make sure the jobs didn't die
    return files


def create_reservation(crispdir, exitneeded=False):
    global resfile
    resfile = op.join(crispdir, 'combine_reservation.sh')
    jobid = os.environ['SLURM_JOB_ID']
    if not op.exists(resfile):
        with open(resfile, 'w') as o:
            o.write("%s" % jobid)
    else:
        exitneeded = True
    time.sleep(random.random()*15)
    with open(resfile, 'r') as o:
        fjobid = o.read().split()[0]
    if not fjobid == jobid or exitneeded is True:
        # just in case two jobs try at nearly the same time
        print('another job has already created reservation for %s' % op.basename(sys.argv[0]))
        exit()
    return jobid


def get_tables(files, pooldir):
    tablefiles = [f for f in fs(op.join(pooldir, 'crisp')) if f.endswith('.txt') and 'all_bedfiles' not in f]
    if not len(tablefiles) == len(files):
        msg = 'for some reason tablefiles != files. jobid=%s' % os.environ['SLURM_JOB_ID']
        print(msg)
        with open(resfile, 'a') as resFile:
            resFile.write("\n%s\n" % msg)
        exit()
    dfs = [remove_multiallelic(thisfile, tablefile, ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)
    
    filename = op.join(pooldir, 'crisp/%s_all_bedfiles_biallelic_snps.txt' % op.basename(pooldir))
    df.to_csv(filename, sep='\t', index=False)
    
    print('combined crisp files to %s' % filename)


def main(pooldir):
    # make sure all of the crisp jobs have finished
    files = checkjobs(pooldir)
    
    # combine table files from output of VariantsToTable
    get_tables(files, pooldir)


if __name__ == "__main__":
    # args
    thisfile, pooldir = sys.argv

    main(pooldir)
