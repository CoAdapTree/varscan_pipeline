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
from balance_queue import getsq
from filter_VariantsToTable import main as remove_multiallelic
from datetime import datetime as dt


def gettimestamp(f):
    return time.ctime(op.getmtime(f))


def getmostrecent(files,remove=False):
    if not isinstance(files, list):
        files = [files]
    if len(files) > 1:
        whichout = files[0]
        dt1 = dt.strptime(gettimestamp(whichout)[4:],"%b %d %H:%M:%S %Y")
        for o in files[1:]:
            dt2 = dt.strptime(gettimestamp(o)[4:],"%b %d %H:%M:%S %Y")
            if dt2 > dt1:
                whichout = o
                dt1 = dt2
        if remove == True:
            [os.remove(f) for f in files if not f == whichout]
        return whichout
    elif len(files) == 1:
        return files[0]
    else:
        # if len(files) == 0
        return None


def checkpids(outs, queue):
    # if any of the other crisp jobs are pending or running, exit
    locals().update({'thisfile': thisfile})
    pids = [q[0] for q in queue]
    jobid = os.environ['SLURM_JOB_ID']
    for out in outs:
        pid = out.split("_")[-1].replace(".out", "")
        if pid in pids and pid != jobid:  # if the job is running, but it's not this job
            print('the following file is still in the queue - exiting %(thisfile)s' % locals(),
                  '\n', '\t%(out)s' % locals())
            exit()  # TODO: should I be exiting?


def check_queue(outs, pooldir):
    # get jobs from the queue, except those that are closing (assumes jobs haven't failed)
    sq = getsq(grepping=['crisp_bedfile', op.basename(pooldir)], states=['R', 'PD'])  
    if len(sq) > 0:
        checkpids(outs, sq)
    # no need for an else statement here, if len(sq) == 0: no need to check the pids

    
def getfiles(samps, shdir, grep):
    found = [sh for sh in fs(shdir) if sh.endswith(".sh") and grep in sh]
    outs = [out for out in fs(shdir) if out.endswith('.out') and grep in out]
    if len(found) != len(samps):
        print('not all shfiles have been created, exiting combine_crisp.py')
        exit()
    files = dict((f, getmostrecent([out for out in outs if op.basename(f).replace(".sh","") in out])) 
                 for f in found)
    if None in files.values():
        print('not all shfiles have been sbatched, exiting combine_crisp.py')
        exit()
    return files


def check_seff(outs):
    for f in outs:
        pid = f.split("_")[-1].replace(".out", "")
        seff, seffcount = '', 0
        while isinstance(seff, list) == False:
            # sometimes slurm sucks
            seff = subprocess.check_output(['seff', pid]).decode('utf-8').split('\n')
            if seffcount == 10:
                print('slurm is screwing something up with seff, exiting')
                exit()
            time.sleep(1)
            seffcount += 1
        state = [x.lower() for x in seff if 'State' in x][0]
        if not 'exit code 0' in state:
            print('job died (%s) for %s' % (state, f)) 
            print('exiting')
            exit()


def checkjobs(pooldir):
    pool = op.basename(pooldir)
    samps = pklload(op.join(op.dirname(pooldir),'poolsamps.pkl'))[pool]
    shdir = op.join(pooldir, 'shfiles/crisp')
    files = getfiles(samps, shdir, 'crisp_bedfile')
#     if not len(files) == len(found):
#         unmade = [f for f in found if f not in files]
#         text = ''
#         for missing in unmade:
#             text = text + "\t%s\n" % missing
#         print('still missing files from these files:\n%(text)s\n%(thisfile)s is exiting\n' % locals())
#         exit()
    check_queue(files.values(), pooldir)  # make sure job isn't in the queue (running or pending)
    check_seff(files.values())  # make sure the jobs didn't die
    return files, shdir


def create_reservation(crispdir, exitneeded=False):
    global resfile, jobid
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
        print('another job has already created reservation for %s' % op.basename(thisfile))
        exit()


def get_tables(files):
    tablefiles = [f for f in fs(op.join(pooldir, 'crisp')) if f.endswith('.txt') and 'all_bedfiles' not in f]
    if not len(tablefiles) == len(files):
        msg = 'for some reason tablefiles != files. jobid=%s' % jobid
        print(msg)
        with open(resfile, 'a') as resFile:
            resFile.write("\n%s\n" % msg)
        exit()
    dfs = [remove_multiallelic(thisfile, tablefile, ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)
    
    filename = op.join(pooldir, 'crisp/%s_all_bedfiles.txt' % op.basename(pooldir))
    df.to_csv(filename, sep='\t', index=False)
    
    print('combined crisp files to %s' % filename)


def main(pooldir):
    # make sure all of the crisp jobs have finished
    files, crispdir = checkjobs(pooldir)
    
    # create reservation so other files don't try and write files.sh, exit() if needed
    create_reservation(crispdir)

    # combine table files from output of VariantsToTable
    get_tables(files)


if __name__ == "__main__":
    # args
    thisfile, pooldir = sys.argv

    main(pooldir)
