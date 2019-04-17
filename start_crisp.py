"""
### purpose
# sbatch crisp cmd if all bamfiles have been created
###

### usage
# python start_crisp.py parentdir pool
###

### fix
# assumes equal sample size across pools
###
"""


import sys, os, time, random, subprocess, balance_queue, shutil
from os import path as op
from coadaptree import makedir, fs, pklload, get_email_info
from balance_queue import getsq


def gettimestamp(f):
    return time.ctime(op.getmtime(f))


def getmostrecent(files, remove=False):
    if not isinstance(files, list):
        files = [files]
    if len(files) > 1:
        whichout = files[0]
        dt1 = dt.strptime(gettimestamp(whichout)[4:], "%b %d %H:%M:%S %Y")
        for o in files[1:]:
            dt2 = dt.strptime(gettimestamp(o)[4:], "%b %d %H:%M:%S %Y")
            if dt2 > dt1:
                whichout = o
                dt1 = dt2
        if remove is True:
            rems = [os.remove(f) for f in files if not f == whichout]
        return whichout
    elif len(files) == 1:
        return files[0]
    else:
        # if len(files) == 0
        return None


def getfiles(samps, shdir, grep):
    if grep == 'crisp_bedfile':
        # numbedfiles = numshfiles = found; so set samps to numbedfiles to bypass len(found) == len(samps)
        pooldir = op.dirname(op.dirname(shdir))
        parentdir = op.dirname(pooldir)
        pool = op.basename(pooldir)
        ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
        samps = fs(op.join(op.dirname(ref), 'bedfiles_%s' % op.basename(ref).split(".fa")[0]))
    found = [sh for sh in fs(shdir) if sh.endswith(".sh") and grep in sh]
    outs = [out for out in fs(shdir) if out.endswith('.out') and grep in out]
    if len(found) != len(samps):
        # this only works for start crisp, not combine crisp since crisp has 15 bedfiles.
        print('not all shfiles have been created, exiting %s' % sys.argv[0])
        exit()
    files = dict((f, getmostrecent([out for out in outs if op.basename(f).replace(".sh", "") in out]))
                 for f in found)
    if None in files.values():
        print('not all shfiles have been sbatched, exiting %s' % sys.argv[0])
        exit()
    return files


def check_seff(outs):
    jobid = os.environ['SLURM_JOB_ID']
    for f in outs:
        pid = f.split("_")[-1].replace(".out", "")
        if not pid == jobid:
            seff, seffcount = '', 0
            while isinstance(seff, list) is False:
                # sometimes slurm sucks
                seff = subprocess.check_output([shutil.which('seff'), pid]).decode('utf-8').split('\n')
                if seffcount == 10:
                    print('slurm is screwing something up with seff, exiting %s' % sys.argv[0])
                    exit()
                time.sleep(1)
                seffcount += 1
            state = [x.lower() for x in seff if 'State' in x][0]
            if 'exit code 0' not in state:
                status = 'died' if 'running' not in state else 'is running'
                print('cannot proceed with %s' % sys.argv[0])
                print('job %s (%s) for %s' % (status, state, f))
                print('exiting %s' % sys.argv[0])
                exit()


def checkpids(outs, queue):
    # if any of the other crisp jobs are pending or running, exit
    pids = [q[0] for q in queue]
    jobid = os.environ['SLURM_JOB_ID']
    for out in outs:
        pid = out.split("_")[-1].replace(".out", "")
        if pid in pids and pid != jobid:  # if the job is running, but it's not this job
            print('the following file is still in the queue - exiting %s' % sys.argv[0],
                  '\n', '\t%(out)s' % locals())
            exit()


def check_queue(outs, pooldir):
    # get jobs from the queue, except those that are closing (assumes jobs haven't failed)
    sq = getsq(grepping=['crisp_bedfile', op.basename(pooldir)], states=['R', 'PD'])
    if len(sq) > 0:
        checkpids(outs, sq)
    # no need for an else statement here, if len(sq) == 0: no need to check the pids


def get_bamfiles(samps, pooldir):
    found = fs(op.join(pooldir, '04_realign'))
    files = dict((samp, f.replace(".bai", ".bam")) for samp in samps for f in found if samp in f and f.endswith('.bai'))
    if not len(files) == len(samps):
        print('len(files) != len(samps)')
        print('files = ', files)
        print('samps = ', samps)
        exit()
    return files


def checkfiles(pooldir):
    # get the list of file names
    pool = op.basename(pooldir)
    samps = pklload(op.join(op.dirname(pooldir), 'poolsamps.pkl'))[pool]
    shdir = op.join(pooldir, 'shfiles/05_indelRealign_shfiles')
    files = getfiles(samps, shdir, 'indelRealign')
    check_queue(files.values(), pooldir)  # make sure job isn't in the queue (running or pending)
    check_seff(files.values())  # make sure the jobs didn't die
    return get_bamfiles(samps, pooldir)
#     samps, files = getfiles()
#     if not len(samps) == len(files):
#         unmade = [samp for samp in samps if samp not in files]
#         text = ''
#         for missing in unmade:
#             text = text + "\t%s\n" % missing
#         print("still missing files from these samps:\n%s\n%s is exiting\n" % (text, thisfile))
#         exit()
#     return list(files.values())


def create_reservation(pooldir, exitneeded=False):
    crispdir = makedir(op.join(pooldir, 'shfiles/crisp'))
    file = op.join(crispdir, '%s_crisp_reservation.sh' % pool)
    jobid = os.environ['SLURM_JOB_ID']
    if not op.exists(file):
        with open(file, 'w') as o:
            o.write("%s" % jobid)
    else:
        exitneeded = True
    time.sleep(random.random()*15)
    with open(file, 'r') as o:
        fjobid = o.read().split()[0]
    if not fjobid == jobid or exitneeded is True:
        # just in case two jobs try at nearly the same time
        print('another job has already created start_crisp.sh for %s' % pool)
        exit()
    return crispdir


def getcmd(bamfiles, bedfile, pool, pooldir, parentdir):
    bams = ' --bam '.join(bamfiles)
    num = bedfile.split("_")[-1].split(".bed")[0]
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    outdir = makedir(op.join(pooldir, 'crisp'))
    outfile = op.join(outdir, '%(pool)s_crisp_bedfile_%(num)s.vcf' % locals())
    poolsize = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
    logfile = outfile.replace(".vcf", ".log")
    convertfile = outfile.replace(".vcf", "_converted.vcf")
    return ('''$CRISP_DIR/CRISP --bam %(bams)s --ref %(ref)s --VCF %(outfile)s \
--poolsize %(poolsize)s --mbq 20 --minc 5 --bed %(bedfile)s > %(logfile)s

touch $SLURM_TMPDIR/bam_file_list.txt # assumes equal pool sizes

$CRISP_DIR/scripts/convert_pooled_vcf.py %(outfile)s $SLURM_TMPDIR/bam_file_list.txt \
%(poolsize)s > %(convertfile)s
''' % locals(),
            num, outfile, convertfile, logfile)


def make_sh(bamfiles, bedfile, crispdir, pool, pooldir):
    cmd, num, outfile, convertfile, logfile = getcmd(bamfiles,
                                                     bedfile,
                                                     pool,
                                                     pooldir,
                                                     op.dirname(pooldir))
    tablefile = convertfile.replace(".vcf", "_table.txt")
    email_text = get_email_info(parentdir, 'crisp')
    text = f'''#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name={pool}-crisp_bedfile_{num}
#SBATCH --time=23:59:00
#SBATCH --mem=15000M
#SBATCH --output={pool}-crisp_bedfile_{num}_%j.out
{email_text}

module load python/2.7.14
# run CRISP (commit 60966e7)
{cmd}
module unload python

# vcf -> table (multiallelic to multiple lines, filtered in combine_crisp.py
module load gatk/4.1.0.0
gatk VariantsToTable --variant {convertfile} -F CHROM -F POS -F REF -F ALT -F AF -F QUAL \
-F DP -F CT -F AC -F VT -F EMstats -F HWEstats -F VF -F VP -F HP -F MQS -F TYPE -F FILTER \
-O {tablefile} --split-multi-allelic
module unload gatk

# gzip outfiles to save space
cd $(dirname {outfile})
gzip {outfile}
gzip {convertfile}
rm {logfile}

# if any other crisp jobs are hanging due to priority, change the account
source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"
python $HOME/pipeline/balance_queue.py crisp

'''
    file = op.join(crispdir, '%(pool)s-crisp_bedfile_%(num)s.sh' % locals())
    with open(file, 'w') as o:
        o.write("%s" % text)
    return file


def sbatch(file):
    os.chdir(op.dirname(file))
    pid = subprocess.check_output([shutil.which('sbatch'), file]).decode('utf-8').replace("\n", "").split()[-1]
    print("sbatched %s" % file)
    return pid


def get_bedfiles(parentdir, pool):
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    beddir = op.join(op.dirname(ref), 'bedfiles_%s' % op.basename(ref).split(".fa")[0])
    return [f for f in fs(beddir) if f.endswith('.bed')]  # TODO: see if I split any other refs by .fasta


def create_sh(bamfiles, crispdir, pool, pooldir):
    bedfiles = get_bedfiles(parentdir, pool)
    pids = []
    for bedfile in bedfiles:
        file = make_sh(bamfiles, bedfile, crispdir, pool, pooldir)
        pids.append(sbatch(file))
    return pids


def create_combine(pids, parentdir, pool):
    pooldir = op.join(parentdir, pool)
    email_text = get_email_info(parentdir, 'crisp')
    if isinstance(pids, list) is not True:
        print('pids is not a list')
        pids = [pids]
    dependencies = '#SBATCH --dependency=afterok:' + ','.join(pids)
    text = f'''#!/bin/bash
#SBATCH --job-name={pool}-combine-crisp
#SBATCH --time=02:59:00
#SBATCH --mem=125000M
#SBATCH --cpus-per-task=1
#SBATCH --output={pool}-combine-crisp_%j.out
{dependencies}
{email_text}


source $HOME/.bashrc
export PYTHONPATH="${{PYTHONPATH}}:$HOME/pipeline"
export SQUEUE_FORMAT="%.8i %.8u %.12a %.68j %.3t %16S %.10L %.5D %.4C %.6b %.7m %N (%r)"

python $HOME/pipeline/combine_crisp.py {pooldir}

'''
    pooldir = op.join(parentdir, pool)
    shdir = op.join(pooldir, 'shfiles/crisp')
    combfile = op.join(shdir, '%s-combine-crisp.sh' % pool)
    with open(combfile, 'w') as o:
        o.write("%s" % text)
    sbatch(combfile)
    print('sbatched combinefile with dependencies: ' + ','.join(pids))


def main(parentdir, pool):
    """Start CRISP if it's appropriate to do so"""

    # check to see if all bam files have been created; if not: exit()
    bamfiles = checkfiles(op.join(parentdir, pool))

    # create reservation so other files don't try and write files.sh, exit() if needed
    crispdir = create_reservation(op.join(parentdir, pool))

    # create .sh file and submit to scheduler
    pids = create_sh(bamfiles.values(), crispdir, pool, op.join(parentdir, pool))

    # create .sh file to combine crisp parallels using jobIDs as dependencies
    create_combine(pids, parentdir, pool)

    # balance queue
    balance_queue.main('balance_queue.py', 'crisp')


if __name__ == "__main__":
    # args
    thisfile, parentdir, pool = sys.argv

    main(parentdir, pool)
