### fix
# getfiles() will run into problems if a job is rescheduled and its old .out file remains in crispdir
###

### purpose
# combine the bedfiles from parallelized crisp runs
###

### usage
# python combine_crisp.py /path/to/pooldir/
###

### imports
import os, sys, time, random, pandas as pd
from os import path as op
from coadaptree import fs, makedir
from balance_queue import checksq
from filter_VariantsToTable import main as remove_multiallelic
###

def getfiles():
    global crispdir
    crispdir = op.join(pooldir,'shfiles/crisp')
    found = [sh for sh in fs(crispdir) if sh.endswith(".sh") and 'crisp_bedfile' in sh]
    outs  = [out for out in fs(crispdir) if out.endswith('.out') and 'crisp_bedfile' in out]
    files = dict((f,out) for f in found for out in outs if op.basename(f).replace(".sh","") in out)
    return (found,files)

def checkpids(files,queue):
    locals().update({'thisfile':thisfile})
    pids  = [x.split()[0] for x in queue]
    jobid = os.environ['SLURM_JOB_ID']
    for f,out in files.items():
        pid = out.split("_")[-1].replace(".out","")
        if pid in pids and pid != jobid: # if the job is running, but it's not this job
            os.system('echo the following file is still in the queue - exiting %(thisfile)s\n\t%(out)s' % locals())
            exit()

def check_queue(files):
    user = os.environ['USER']
    # get jobs from the queue, except those that are closing
    SQ = os.popen('''squeue -u %(user)s | grep "crisp_bedfile" | grep -v "CG"''' % locals()).read().split("\n")
    SQ = [s for s in SQ if not s == '']
    if len(SQ) > 0:
        checksq(SQ,'combine_crisp')
        checkpids(files,SQ)

def checkjobs():
    locals().update({'thisfile':thisfile})
    found,files = getfiles()
    if not len(files) == len(found):
        unmade = [f for f in found if not f in files]
        text = ''
        for missing in unmade:
            text = text + "\t%s\n" % missing
        os.system('echo still missing files from these files:\n%(text)s\n%(thisfile)s is exiting\n' % locals())
        exit()
    check_queue(files) # make sure job isn't in the queue
    return files
        
def create_reservation(exitneeded=False):
    global resfile, jobid
    resfile = op.join(crispdir,'combine_reservation.sh')
    jobid = os.environ['SLURM_JOB_ID']
    if not op.exists(resfile):
        with open(resfile,'w') as o:
            o.write("%s" % jobid)
    else:
        exitneeded = True
    time.sleep(random.random()*15)
    with open(resfile,'r') as o:
        fjobid = o.read().split()[0]
    if not fjobid == jobid or exitneeded == True:
        # just in case two jobs try at nearly the same time
        os.system('echo another job has already created reservation for %s' % op.basename(thisfile) )
        exit()
        
def get_tables(files):
    tablefiles = [f for f in fs(op.join(pooldir,'crisp')) if f.endswith('.txt') and 'all_bedfiles' not in f]
    if not len(tablefiles) == len(files):
        msg = 'for some reason tablefiles != files. jobid=%s' % jobid
        os.system('echo %s' % msg)
        os.system('echo -e "\n%s" >> %s' % (msg,resfile)) # it's ok to write a msg, this will be the only time to try combine
        exit()
    dfs = [remove_multiallelic(thisfile,tablefile,ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)
    
    filename = op.join(pooldir,'crisp/%s_all_bedfiles.txt' % op.basename(pooldir))
    df.to_csv(filename,sep='\t',index=False)
    
    os.system('echo combined crisp files to %s' % filename)

def main():
    # make sure all of the crisp jobs have finished
    files = checkjobs()
    
    # create reservation so other files don't try and write files.sh, exit() if needed
    create_reservation()
    
    # combine table files from output of VariantsToTable
    get_tables(files)
    


if __name__ == "__main__":
    # args
    thisfile, pooldir = sys.argv

    main()
    