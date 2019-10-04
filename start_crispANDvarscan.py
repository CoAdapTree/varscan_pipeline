"""Create and sbatch crisp and varscan command files if all realigned bamfiles have been created.

If a single sample per pool_name:
    - set min freq to 0
If multiple samps per pool_name:
    - set min freq to 1/(ploidy_per_samp * nsamps)

# usage
# python start_crispANDvarscan.py parentdir pool
#

# fix
# assumes equal sample size across pools
#
"""


import sys, os, time, random, subprocess, shutil
from os import path as op
from datetime import datetime as dt
from coadaptree import makedir, fs, pklload, get_email_info
from balance_queue import getsq


def gettimestamp(f):
    """Get last time modified."""
    return time.ctime(op.getmtime(f))


def getmostrecent(files):
    """Determine the most recent file from a list of files."""
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
        return whichout
    elif len(files) == 1:
        return files[0]
    else:
        # if len(files) == 0
        return None


def getfiles(samps, shdir, grep):
    """Determine if all realign bam jobs have been created and sbatched.

    Positional arguments:
    samps - list of sample names (length is number of expected shfiles)
    shdir - directory where .sh and .out files are
    grep - program name - keyword used to find correct files

    Returns:
    files - dictionary where key = sh file, val = most recent outfile
    """
    found = [sh for sh in fs(shdir) if sh.endswith(".sh") and grep in sh]
    outs = [out for out in fs(shdir) if out.endswith('.out') and grep in out]
    if len(found) != len(samps):
        print('not all shfiles have been created, exiting %s' % sys.argv[0])
        exit()
    files = dict((f, getmostrecent([out for out in outs if op.basename(f).replace(".sh", "") in out]))
                 for f in found)
    if None in files.values():
        print('not all shfiles have been sbatched, exiting %s' % sys.argv[0])
        exit()
    return files


def check_seff(outs):
    """Execute slurm seff command on each outfile's slurm_job_id to ensure it ran without error.
    Exit otherwise.
    """
    print('checking seff')
    jobid = os.environ['SLURM_JOB_ID']
    for i,f in enumerate(outs):
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
        if (i+1) % 10 == 0:
            print('\t%s/%s' % (i+1, len(outs)))


def checkpids(outs, queue):
    """If any of the other crisp jobs are pending or running, exit."""
    print('checking pids')
    pids = [q[0] for q in queue]
    jobid = os.environ['SLURM_JOB_ID']
    for out in outs:
        pid = out.split("_")[-1].replace(".out", "")
        if pid in pids and pid != jobid:  # if the job is running, but it's not this job
            print('the following file is still in the queue - exiting %s' % sys.argv[0],
                  '\n', '\t%(out)s' % locals())
            exit()


def check_queue(outs, pooldir):
    """Get jobs from the queue, except those that are closing (assumes jobs haven't failed)."""
    print('checking queue')
    sq = getsq(grepping=['crisp_bedfile', op.basename(pooldir)])
    if len(sq) > 0:
        checkpids(outs, sq)
    # no need for an else statement here, if len(sq) == 0: no need to check the pids


def get_bamfiles(samps, pooldir):
    """Using a list of sample names, find the realigned bamfiels.

    Return:
    files - dictionary with key = samp_name, val = /path/to/bamfile
    """
    print('getting bamfiles')
    found = fs(op.join(pooldir, '04_realign'))
    files = dict((samp, f.replace(".bai", ".bam")) for samp in samps for f in found if samp in f and f.endswith('.bai'))
    if not len(files) == len(samps):
        print('len(files) != len(samps)')
        print('files = ', files)
        print('samps = ', samps)
        exit()
    return files


def checkfiles(pooldir):
    """Call get_bamfiles."""
    # get the list of file names
    print('checking files')
    pool = op.basename(pooldir)
    samps = pklload(op.join(op.dirname(pooldir), 'poolsamps.pkl'))[pool]
    shdir = op.join(pooldir, 'shfiles/05_indelRealign_shfiles')
    files = getfiles(samps, shdir, 'indelRealign')
    check_queue(files.values(), pooldir)  # make sure job isn't in the queue (running or pending)
    check_seff(files.values())  # make sure the jobs didn't die
    return get_bamfiles(samps, pooldir)


def create_reservation(pooldir, exitneeded=False):
    """Create a file so that other realign jobs can't start crisp and varscan too."""
    print('creating reservation')
    shdir = makedir(op.join(pooldir, 'shfiles/crispANDvarscan'))
    file = op.join(shdir, '%s_crispANDvarscan_reservation.sh' % pool)
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
        print('\tanother job has already created crispANDvarscan_reservation.sh for %s' % pool)
        exit()
    return shdir


def get_prereqs(bedfile, parentdir, pool, program):
    """Get object names."""
    num = bedfile.split("_")[-1].split(".bed")[0]
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
#     outdir = makedir(op.join(pooldir, program))
#     vcf = op.join(outdir, f'{pool}_{program}_bedfile_{num}.vcf')
    vcf = f'$SLURM_TMPDIR/{pool}_{program}_bedfile_{num}.vcf'
    return (num, ref, vcf)


def get_small_bam_cmds(bamfiles, bednum, bedfile):
    """Get samtools commands to reduce a bamfile to intervals in the bedfile."""
    smallbams = []
    cmds = '''module load java\nmodule load samtools/1.9\n'''
    for bam in bamfiles:
        pool = op.basename(bam).split("_realigned")[0]
        smallbam = f'$SLURM_TMPDIR/{pool}_realigned_{bednum}.bam'
        cmd = f'''samtools view -b -L {bedfile} {bam} > {smallbam}\n'''
        cmds = cmds + cmd
        smallbams.append(smallbam)
    return (smallbams, cmds)


def get_crisp_cmd(bamfiles, bedfile, pool, parentdir, ref, vcf, bednum):
    """Create command to call crisp."""
    smallbams, smallcmds = get_small_bam_cmds(bamfiles, bednum, bedfile)
    bams = ' --bam '.join(smallbams)
    poolsize = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
    logfile = vcf.replace(".vcf", ".log")
    convertfile = vcf.replace(".vcf", "_converted.vcf")
    cmds = smallcmds + f'''module load python/2.7.14
$CRISP_DIR/CRISP --bam {bams} --ref {ref} --VCF {vcf} \
--poolsize {poolsize} --mbq 20 --minc 5 --bed {bedfile} > {logfile}

touch $SLURM_TMPDIR/bam_file_list.txt # assumes equal pool sizes

$CRISP_DIR/scripts/convert_pooled_vcf.py {vcf} $SLURM_TMPDIR/bam_file_list.txt \
{poolsize} > {convertfile}
module unload python
'''
    return cmds, convertfile, logfile


def check_filtering_options(pooldir, vcf, finalvcf, filtering_cmd=''):
    """
    See if the user specified any additional filtering during 00_start.py.

    Assumes varscan as program.
    """
    nextvcf = vcf # will be overwritten if filtered; vcf is in SLURM_TMPDIR anway
    parentdir = op.dirname(pooldir)
    # check for paralog filtering
    parafile = op.join(parentdir, 'paralog_snps.pkl')
    if op.exists(parafile):
        paralogs = pklload(parafile)
        nextvcf = vcf.replace('.vcf', '_paralog-filtered.vcf')
        paralogsnps = finalvcf.replace(".vcf", "_paralogs.vcf.gz")
        filtering_cmd = filtering_cmd + f'''# Remove paralog sites
vcftools --vcf {vcf} --out {nextvcf} --exclude-positions {paralogs} --recode --recode-INFO-all
# put paralog SNPs into file for later use
vcftools --vcf {vcf} --positions {paralogs} --stdout --recode --recode-INFO-all | gzip -c > {paralogsnps}
'''
        vcf = nextvcf
    # check for repeat masking
    repeatfile = op.join(parentdir, 'repeat_regions.pkl')
    if op.exists(repeatfile):
        repeats = pklload(repeatfile)
        nextvcf = vcf.replace('.vcf', '_repeat-filtered.vcf')
        filtering_cmd = filtering_cmd + f'''# Remove repeat-masked regions
vcftools --vcf {vcf} --out {nextvcf} --exclude-bed {repeats} --recode --recode-INFO-all
'''
    # create final fintering_cmd
    if filtering_cmd == '':
        # if no filters were applied, move the varscan output from tmp to parentsubdir
        filtering_cmd = f'''
# move file
mv {nextvcf} {finalvcf}
'''
    else:
        filtering_cmd = f'''
# begin filtering vcf files
module load vcftools/0.1.14
{filtering_cmd}
module unload vcftools
# move file
mv {nextvcf} {finalvcf}
'''
    return filtering_cmd

def get_varscan_cmd(bamfiles, bedfile, bednum, vcf, ref, pooldir, program):
    """Create command to call varscan."""
    smallbams, smallcmds = get_small_bam_cmds(bamfiles, bednum, bedfile)
    smallbams = ' '.join(smallbams)
    ploidy = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
    # if single-sample then set minfreq to 0, else use min possible allele freq
    minfreq = 1/(ploidy*len(bamfiles)) if len(bamfiles) > 1 else 0
    cmd = f'''samtools mpileup -B -f {ref} {smallbams} | java -Xmx15g -jar \
$VARSCAN_DIR/VarScan.v2.4.3.jar mpileup2cns --min-coverage 8 --p-value 0.05 \
--min-var-freq {minfreq} --strand-filter 1 --min-freq-for-hom 0.80 \
--min-avg-qual 20 --output-vcf 1 > {vcf}
module unload samtools
'''
    # final vcf
    outdir = makedir(op.join(pooldir, program))
    finalvcf = op.join(outdir, op.basename(vcf))
    filtering_cmd = check_filtering_options(pooldir, vcf, finalvcf)
    cmds = smallcmds + cmd + filtering_cmd
    return (cmds, finalvcf)


def make_sh(bamfiles, bedfile, shdir, pool, pooldir, program, parentdir):
    """Create sh file for varscan or crisp command."""
    num, ref, vcf = get_prereqs(bedfile, parentdir, pool, program)
    if program == 'crisp':
        cmd, finalvcf, logfile = get_crisp_cmd(bamfiles,
                                               bedfile,
                                               pool,
                                               parentdir,
                                               ref,
                                               vcf,
                                               num)
        # force gzip the vcf (overwrite previous .gz), rm logfile
        second_cmd = f'''bgzip -f {vcf}
rm {logfile}
'''
        mem = "9000M"
        time = '3-00:00:00'
        fields = '''-F DP -F CT -F AC -F VT -F EMstats -F HWEstats -F VF -F VP \
-F HP -F MQS -GF GT -GF GQ -GF DP'''
    else:
        cmd, finalvcf = get_varscan_cmd(bamfiles, bedfile, num,
                                        vcf, ref, pooldir, program)
        second_cmd = ''''''
        mem = "2000M"
        time = '7-00:00:00'
        fields = '''-F ADP -F WT -F HET -F HOM -F NC -GF GT -GF GQ -GF SDP -GF DP \
-GF FREQ -GF PVAL -GF AD'''

    tablefile = finalvcf.replace(".vcf", "_table.txt")
    bash_variables = op.join(parentdir, 'bash_variables')
    text = f'''#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --job-name={pool}-{program}_bedfile_{num}
#SBATCH --time={time}
#SBATCH --mem={mem}
#SBATCH --output={pool}-{program}_bedfile_{num}_%j.out

# run CRISP (commit 60966e7) or VarScan (v.2.4.2)
{cmd}

# vcf -> table (multiallelic to multiple lines, filtered in combine_crispORvarscan.py
module load gatk/4.1.0.0
gatk VariantsToTable --variant {finalvcf} -F CHROM -F POS -F REF -F ALT -F AF -F QUAL \
-F TYPE -F FILTER {fields} -O {tablefile} --split-multi-allelic
module unload gatk

# gzip outfiles to save space
module load nixpkgs/16.09  gcc/7.3.0 htslib/1.9
cd $(dirname {finalvcf})
bgzip -f {finalvcf}
{second_cmd}

# if any other varscan jobs are hanging due to priority, change the account
source {bash_variables}
python $HOME/pipeline/balance_queue.py {program} {parentdir}

'''
    file = op.join(shdir, f'{pool}-{program}_bedfile_{num}.sh')
    with open(file, 'w') as o:
        o.write("%s" % text)
    return file


def sbatch(file):
    """Sbatch file."""
    os.chdir(op.dirname(file))
    pid = subprocess.check_output([shutil.which('sbatch'), file]).decode('utf-8').replace("\n", "").split()[-1]
    print("sbatched %s" % file)
    time.sleep(0.25)
    return pid


def get_bedfiles(parentdir, pool):
    """Get a list of paths to all of the bed files for ref.fa."""
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    beddir = op.join(op.dirname(ref), 'bedfiles_%s' % op.basename(ref).split(".fa")[0])
    return [f for f in fs(beddir) if f.endswith('.bed')]


def create_sh(bamfiles, shdir, pool, pooldir, program, parentdir):
    """Create and sbatch shfiles, record pid to use as dependency for combine job."""
    bedfiles = get_bedfiles(parentdir, pool)
    pids = []
    for bedfile in bedfiles:
        file = make_sh(bamfiles, bedfile, shdir, pool, pooldir, program, parentdir)
        pids.append(sbatch(file))
    return pids


def create_combine(pids, parentdir, pool, program, shdir):
    """Create command file to combine crisp or varscan jobs once they're finished.

    Positional arguments:
    pids = list of slurm job id dependencies (the jobs that need to finish first)
    ...
    """
    pooldir = op.join(parentdir, pool)
    email_text = get_email_info(parentdir, 'final')
    dependencies = '#SBATCH --dependency=afterok:' + ','.join(pids)
    bash_variables = op.join(parentdir, 'bash_variables')
    text = f'''#!/bin/bash
#SBATCH --job-name={pool}-combine-{program}
#SBATCH --time=12:00:00
#SBATCH --mem=20000M
#SBATCH --cpus-per-task=1
#SBATCH --output={pool}-combine-{program}_%j.out
{dependencies}
{email_text}


source {bash_variables}

python $HOME/pipeline/combine_crispORvarscan.py {pooldir} {program} {pool}

'''
    combfile = op.join(shdir, f'{pool}-combine-{program}.sh')
    with open(combfile, 'w') as o:
        o.write("%s" % text)
    sbatch(combfile)
    print(f'sbatched {program} combinefile with dependencies: ' + ','.join(pids))


def main(parentdir, pool):
    """Start <program> if it's appropriate to do so."""

    # check to see if all bam files have been created; if not: exit()
    bamfiles = checkfiles(op.join(parentdir, pool))

    # create reservation so other files don't try and write files.sh, exit() if needed
    shdir = create_reservation(op.join(parentdir, pool))

    # create .sh files
    #for program in ['crisp', 'varscan']:  # I'll be deprecating crisp soon
    for program in ['varscan']:
        print('starting %s commands' % program)
        # create .sh file and submit to scheduler
        pids = create_sh(bamfiles.values(),
                         shdir,
                         pool,
                         op.join(parentdir, pool),
                         program,
                         parentdir)

        # create .sh file to combine crisp parallels using jobIDs as dependencies
        create_combine(pids, parentdir, pool, program, shdir)


if __name__ == "__main__":
    # args
    thisfile, parentdir, pool = sys.argv

    main(parentdir, pool)
