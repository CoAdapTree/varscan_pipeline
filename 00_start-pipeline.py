"""
start the pipeline.

### usage
# 00_start-pipeline.py -p PARENTDIR
#                     [-e EMAIL [-n EMAIL_OPTIONS ...]]
#                     [-maf MAF]
#                     [--rm_paralogs]
#                     [--rm_repeats]
#                     [--translate]
###

### assumes
# that samples duplicated in the 'sample_name' column have the same rglb, rgsm, and rgpl read groups
# that ploidy is the same (ie uses the exact same pop members) across resquencing of ...
#     individual pools (all resequencing must have same ploidy for any given sample)
###

### TODO
# use @functools.wraps(function) for handling [--rm_repeats/paralogs, --translate] to simplify code
"""

import os, sys, distutils.spawn, subprocess, shutil, argparse, pandas as pd
import balance_queue, create_bedfiles
from os import path as op
from subprocess import PIPE
from subprocess import Popen
from collections import OrderedDict
from coadaptree import fs, pkldump, uni, makedir, askforinput, Bcolors, luni

def get_rgid(r1):
    """If RGID is blank, print out the output that the pipeline would otherwise use."""
    p1 = Popen(['zcat',r1], stdout=PIPE)
    p2 = Popen(['head', '-n1'], stdin=p1.stdout, stdout=PIPE)
    return '_'.join(p2.communicate()[0].split()[0].decode('utf-8').split('\n')[0].split(":")[:4])


def create_sh(pooldirs, poolref, parentdir):
    """Run 01_trim-fastq.py to sbatch trimming jobs, then balance queue.

    Positional arguments:
    pooldirs - a list of subdirectories in parentdir for gsroups of pools
    poolref - dictionary with key = pool, val = /path/to/ref
    """
    # create sh files
    print(Bcolors.BOLD + '\nwriting sh files' + Bcolors.ENDC)
    for pooldir in pooldirs:
        pool = op.basename(pooldir)
        print(Bcolors.BOLD + '\npool = %s' % pool + Bcolors.ENDC)
        ref = poolref[pool]
        print('\tsending pooldir and ref to 01_trim-fastq.py')
        subprocess.call([shutil.which('python'),
                         op.join(os.environ['HOME'], 'pipeline/01_trim-fastq.py'),
                         pooldir,
                         ref])
    print("\n")
    balance_queue.main('balance_queue.py', 'trim', parentdir)


def get_datafiles(parentdir, f2pool, data):
    """Get list of files from datatable, make sure they exist in parentdir.
    Create symlinks in /parentdir/<pool_name>/.

    Positional arguments:
    parentdir - directory with datatable.txt and (symlinks to) fastq data
    f2pool - dictionary with key = file.fastq, val = pool_name
    data - datatable.txt with info for pipeline
    """
    print(Bcolors.BOLD + '\nchecking for existance of fastq files in datatable.txt' + Bcolors.ENDC)
    files = [f for f in fs(parentdir) if 'fastq' in f and 'md5' not in f]
    datafiles = data['file_name_r1'].tolist()
    for x in data['file_name_r2'].tolist():
        datafiles.append(x)
    if len(files) > len(datafiles):
        desc = 'more'
    if len(files) < len(datafiles):
        desc = 'less'
    try:
        print(Bcolors.WARNING +
              'WARN: there are %s fastq files in %s than in datatable.txt' % (desc, parentdir) +
              Bcolors.ENDC)
        print(Bcolors.BOLD + 'Here are the files in %s' % parentdir + Bcolors.ENDC)
        for x in files:
            print(op.basename(x))
        print(Bcolors.BOLD + 'Here are the files in datatable.txt' + Bcolors.ENDC)
        for x in datafiles:
            print(x)
        askforinput(newline='')

    except NameError:
        pass

    # create symlinks in pooldirs for visualization
    for f in datafiles:
        src = op.join(parentdir, f)
        if not op.exists(src):
            # make sure file in datatable exists
            print("could not find %s in %s\nmake sure file_name in datatable is its basename" % (f, parentdir))
            print("(symlinks in parentdir to fastq files in other dirs works fine, and is the intentional use)")
            sys.exit(1)
        pooldir = op.join(parentdir, f2pool[f])
        dst = op.join(pooldir, f)
        if not op.exists(dst):
            # easy to visualize in cmdline if script is finding correct group of files by ls-ing pooldir
            os.symlink(src, dst)

    # print out RGID if RGID is none


def make_pooldirs(data, parentdir):
    """Create subdirectories of parentdir.

    Positional arguments:
    data - datatable.txt with info for pipeline
    parentdir - directory with datatable.txt and (symlinks to) fastq data
    """
    # make pool dirs
    print(Bcolors.BOLD + "\nmaking pool dirs" + Bcolors.ENDC)
    pools = uni(data['pool_name'].tolist())
    pooldirs = []
    for p in pools:
        pooldir = op.join(parentdir, p)
        if op.exists(pooldir):
            text = "\tWARN: The pooldir already exists, this WILL overwrite and/or delete previous data: %s" % pooldir
            print(Bcolors.WARNING + text + Bcolors.ENDC)
            askforinput(tab='\t', newline='')
            # first unlink fastq files
            for f in fs(pooldir):
                if f.endswith('.gz'):
                    os.unlink(f)
            # then just delete the directory
            shutil.rmtree(pooldir)
        pooldirs.append(makedir(pooldir))
    return pooldirs


def create_all_bedfiles(poolref, numpools):
    """For each unique ref.fa in datatable.txt, create bedfiles for varscan.

    Positional arguments:
    poolref - dictionary with key = pool, val = /path/to/ref
    """
    # create bedfiles for varscan
    print(Bcolors.BOLD + "\ncreating bedfiles" + Bcolors.ENDC)
    for ref in uni(poolref.values()):
        create_bedfiles.main(ref, numpools)

def choose_file(files, pool, purpose, keep=None):
    """Choose which repeat/paralog file to use if multiple exist."""
    print(Bcolors.BOLD + '\t\tWhich file would you like to use for pool (%s) to %s?' % (pool,purpose) + Bcolors.ENDC)
    nums = []
    for i,f in enumerate(files):
        print('\t\t\t%s %s' % (i,op.basename(f)))
        nums.append(i)
    while True:
        inp = int(input(Bcolors.WARNING + "\t\tINPUT NEEDED: Choose file by number: " + Bcolors.ENDC).lower())
        if inp in nums:
            keep = files[inp]
            break
        else:
            print(Bcolors.FAIL + "Please respond with a number from above." + Bcolors.ENDC)
    # make sure they've chosen at least one account
    while keep is None:
        print(Bcolors.FAIL + "FAIL: You need to specify at least file. Revisiting files..." + Bcolors.ENDC)
        keep = choose_file(files, pool, purpose, keep=None)
    return keep


def get_parafile(parentdir, pool):
    """Obtain file containing paralog SNPs to be removed from final SNPs."""
    parafiles = [f for f in fs(parentdir) if f.endswith('_paralog_snps.txt')]
    if len(parafiles) > 1:
        parafile = choose_file(parafiles, pool, 'remove paralogs')
    elif len(parafiles) == 0:
        parafile = None
    elif len(parafiles) == 1:
        parafile = parafiles[0]
    return parafile

def check_ref_assumptions(samp, ref):
    """Make sure pipeline assumptions about ref.fasta are met."""
    if not op.exists(ref):
        text = 'FAIL: ref for %s does not exist in path: %s' % (samp, ref)
        print(Bcolors.FAIL + text + Bcolors.ENDC)
        print('exiting 00_start-pipeline.py')
        exit()
    needed = []
    for suffix in ['.dict', '.amb', '.ann', '.bwt', '.fai', '.pac', '.sa']:
        refext = ref + suffix if suffix != '.dict' else ref.split('.fa')[0] + suffix
        if not op.exists(refext):
            needed.append(refext)
    if len(needed) > 0:
        print(Bcolors.FAIL +
              'FAIL: the following extensions of the reference are needed to continue, \
please create these files' +
              Bcolors.ENDC)
        for n in needed:
            print(Bcolors.FAIL + n + Bcolors.ENDC)
        print('exiting')
        exit()
    return ref


def handle_repeats(repeats, pool2repeatsfile, ref, data, pool):
    """If --rm_repeats flag was used, determine which repeats file to use."""
    repeatsfile = None
    if repeats is True:
        if pool not in pool2repeatsfile:
            # if more than one pool ask user
            repeatfile = ref.split(".fa")[0] + '_repeats.txt'
            if luni(data.loc[data['ref']==ref, 'pool_name']) > 1:
                # if there is more than one pool with the same ref, confirm repeats is for this pool
                inp = askforinput(tab='\t',
                                  msg='Would you like to use this file to remove repeats from the %s pool?: \n\t\t%s'
                                  % (pool, repeatfile),
                                  newline='')
            else:
                # otherwise there is only one pool with this ref - use the repeats
                inp = 'yes'
            if inp == 'yes':
                if not op.exists(repeatfile):
                    text = 'FAIL: You have indicated that you would like repeat regions removed. \n'
                    text = text + 'FAIL: But the file expected by the pipeline (the one you chose above) \n'
                    text = text + 'FAIL: does not exist: %s' % repeatfile
                    print(Bcolors.FAIL + text + Bcolors.ENDC)
                    exit()
                repeatsfile = repeatfile
    return repeatsfile


def handle_translate(translate, pool2translate, ref, data, pool):
    """If --translate flag was used, determine which file to use for tranlsation."""
    transfile = None
    if translate is True:
        if pool not in pool2translate:
            orderfile = ref.split(".fa")[0] + '.order'
            if luni(data.loc[data['ref']==ref, 'pool_name']) > 1:
                # if there is more than one pool with the same ref, confirm repeats is for this pool
                inp = askforinput(tab='\t',
                                  msg='Would you like to use this file to translate stitched positions for the %s pool?:\
\n\t\t%s' % (pool, orderfile),
                                  newline='')
            else:
                # otherwise there is only one pool with this ref - use the repeats
                inp = 'yes'
            if inp == 'yes':
                if not op.exists(orderfile):
                    text = 'FAIL: You have indicated that you would like stitched regions translated. \n'
                    text = text + 'FAIL: But the file expected by the pipeline (the one you chose above) \n'
                    text = text + 'FAIL: does not exist: %s' % orderfile
                    text = text + '\nexiting 00_start-pipeline.py'
                    print(Bcolors.FAIL + text + Bcolors.ENDC)
                    exit()
                # if the user says they want this pool translated, store orderfile
                transfile = orderfile
    return transfile


def handle_paralogs(paralogs, pool2paralogfile, data, pool, parentdir):
    """If --rm_paralogs was used, determine which file to use."""
    parafile = None
    if paralogs is True:
        if pool not in pool2paralogfile:
            if luni(data['pool_name']) > 1:
                # if more than one pool, see if --rm_paralogs applies to this pool
                inp = askforinput(tab='\t',
                                  msg='Would you like to remove paralogs from the pool: %s?' % pool,
                                  newline='')
                if inp == 'yes':
                    parafile = get_parafile(parentdir, pool)
            else:
                # try to assign paralog file manually, ask if necessary
                parafile = get_parafile(parentdir, pool)
    return parafile


def read_datatable(parentdir):
    """Read in datatable.txt."""
    # look for datatable.txt in parentdir
    datatable = op.join(parentdir, 'datatable.txt')
    if not op.exists(datatable):
        print(Bcolors.FAIL + '''FAIL: the datatable is not in the necessary path: %s
FAIL: exiting 00_start-pipeline.py''' % datatable + Bcolors.ENDC)
        sys.exit(3)
    data = pd.read_csv(datatable, sep='\t')
    return data


def handle_dict_fails(pool2repeatsfile, pool2translate, pool2paralogfile, repeats, translate, paralogs, data, parentdir):
    flagexit = False
    for dic,flag,word in zip([pool2repeatsfile, pool2translate, pool2paralogfile],
                             [repeats, translate, paralogs],
                             ['remove repeats', 'translate stitched positions', 'remove paralogs']):
        # if a flag was specified but none of the pools were selected:
        if flag is True and sum([1 for v in dic.values() if v is not None])==0:
            flagexit = True
            text = 'FAIL: You have indicated that you would like to %s from final SNPs.\n' % word
            text = text + 'FAIL: But the user has not specified at least one pool to %s. \n' % word
            text = text + 'FAIL: You need to respond "yes" to at least one of the prompts above \n'
            text = text + 'FAIL: for assigning a file to a pool - i.e., to use the \n'
            text = text + 'FAIL: %s flag, you must apply it to at least one pool. \n' % word
            if 'repeats' in word:
                text = text + 'FAIL: The file containing repeat regions should be one of the following:\n'
                for ref in uni(data['ref']):
                    repeatfile = ref.split(".fa")[0] + '_repeats.txt'
                    text = text + "\t %s \n" % repeatfile
            elif 'stitched' in word:
                text = text + 'FAIL: The file to translate stitched to unstitched positions should \n'
                text = text + 'FAIL: be one of the following:\n'
                for ref in uni(data['ref']):
                    orderfile = ref.split(".fa")[0] + '.order'
                    text = text + "\t %s \n" % orderfile
            elif 'paralogs' in word:
                text = text + 'FAIL: The file(s) to remove paralogs must be in %s \n' % parentdir
                text = text + 'FAIL: and end with "_paralog_snps.txt".'
            print(Bcolors.FAIL + text + Bcolors.ENDC)
    if flagexit is True:
        exit()


def handle_rg_fails(failing, warning, parentdir, data):
    if len(failing) > 0:
        print(Bcolors.FAIL + 'FAIL: The following samples have blank RG info.' + Bcolors.ENDC)
        for fail in failing:
            print(Bcolors.FAIL + "FAIL: %s" % fail + Bcolors.ENDC)
        print('exiting 00_start-pipeline.py')
        exit()
    if len(warning) > 0:
        outputs = []
        for row in data.index:
            samp = data.loc[row, 'sample_name']
            if samp in warning:
                r1 = op.join(parentdir, data.loc[row, 'file_name_r1'])
                outputs.append("\t\t%s\t%s" % (samp, get_rgid(r1)))
        print(Bcolors.WARNING + '\n\n\tWARN: at least one of the samples has a blank RGID in the datatable.\n' +
              '\tWARN: If RGPU is also blank, the pipeline will assign RGPU as: $RGID.$RGLB\n' +
              '\tWARN: The pipeline will automatically assign the following RGIDs.\n' +
              '\n\t\tsample_name\tassigned_RGID' +
              Bcolors.ENDC)
        for output in outputs:
            print(Bcolors.WARNING + output + Bcolors.ENDC)
        askforinput(tab='\t', newline='')


def parse_datatable(data, parentdir, translate, repeats, paralogs):
    """
    Checks some assumptions of datatable.txt, create files and dirs for downstream.

    translate, repeats, and paralogs are boolean.
    parentdir is a path.
    """
    print(Bcolors.BOLD + '\nReading datatable, getting fastq info' + Bcolors.ENDC)
    
    # inititate dictionaries for downstream pipeline
    rginfo = {}  # key=samp vals=rginfo
    samp2pool = {}  # key=samp val=pool
    poolref = {}  # key=pool val=ref.fa
    ploidy = {}  # key=pool val=dict(key=sample: val=sample_ploidy)
    poolsamps = {}  # key=pool val=sampnames
    f2samp = {}  # key=f val=samp
    f2pool = {}  # key=f val=pool
    adaptors = OrderedDict()  # key=samp val={'r1','r2'} val=adaptor
    warning = []  # whether to print out warning about optional RG info
    failing = []  # whether to print out failing about required RG info
    pool2paralogfile = {}  # if --rm_paralogs flagged, store file based on pool
    pool2repeatsfile = {}  # if --rm_repeats flagged, store file based on pool
    pool2translate = {}  # if --translate flagged, store file based on pool

    # make sure there are no blanks where there shouldn't be
    badcols = []
    for column in data.columns:
        if column not in ['rgid', 'rgpu', 'adaptor_1', 'adaptor_2']:
            if data[column].isnull().sum() > 0:
                badcols.append(column)
    if len(badcols) > 0:
        print(Bcolors.FAIL + "\tFAIL: Some rows in datable.txt have blank entries in the following columns: " + Bcolors.ENDC)
        for col in badcols:
            print(Bcolors.FAIL + "\tFAIL: %s" % col + Bcolors.ENDC)
        print('exiting 00_start-pipeline.py')
        exit()

    # make sure specific words are not in a pool name
    badnames = []
    for pool in uni(data['pool_name']):
        for keyword in ['SNP', 'REPEAT', 'PARALOG']:
            if keyword in pool:
                badnames.append((pool, keyword))
    if len(badnames) > 0:
        print(Bcolors.FAIL + "\tFAIL: Some pool names have characters that could cause errors downstream." + Bcolors.ENDC)
        print(Bcolors.FAIL + "\tFAIL: Remove the bad characters from pool_names to continue." + Bcolors.ENDC)
        for pool,keyword in badnames:
            print(Bcolors.FAIL + "\tFAIL: Remove '%s' from pool_name '%s'." % (keyword, pool))
        print('exiting 00_start-pipeline.py')
        exit()

    # iterate through datatable
    for row in data.index:
        # get variables
        samp = data.loc[row, 'sample_name']
        adaptors[samp] = {'r1': data.loc[row, 'adaptor_1'],
                          'r2': data.loc[row, 'adaptor_2']}
        pool = data.loc[row, 'pool_name']
        pooldir = op.join(parentdir, pool)
        print('\t{}\tsamp = {}\tpool = {}'.format(row, samp, pool))
        if pool not in poolsamps:
            poolsamps[pool] = []
        if samp not in poolsamps[pool]:
            poolsamps[pool].append(samp)
        if samp in samp2pool:
            if samp2pool[samp] != pool:
                print(Bcolors.FAIL + 'FAIL: there are duplicate sample names with \
different pool assignments: %s' % samp + Bcolors.ENDC)
                print('exiting')
                exit()
        samp2pool[samp] = pool

        # get ploidy info
        if pool not in ploidy:
            ploidy[pool] = {}
        if samp in ploidy[pool].keys():
            if ploidy[pool][samp] != int(data.loc[row, 'ploidy']):
                text = "FAIL: the ploidy values for sample_name '%s' are not the same" % samp
                print(Bcolors.FAIL + text + Bcolors.ENDC)
                exit()
        ploidy[pool][samp] = int(data.loc[row, 'ploidy'])

        # get ref.fasta info
        ref = data.loc[row, 'ref']
        if pool in poolref:
            # make sure each row for a pool specifies the same reference.fa
            if poolref[pool] != ref:
                text = "FAIL: Ref genome for samples in %s pool seem to have different paths in datatable" % pool
                print(Bcolors.FAIL + text + Bcolors.ENDC)
                print('exiting 00_start-pipeline.py')
                exit()
        else:
            # check assumptions about ref
            poolref[pool] = check_ref_assumptions(samp, ref)

        # hangle RG info
        rginfo[samp] = {}
        # required RG info
        for col in ['rglb', 'rgpl', 'rgsm']:  # rg info columns
            if not data.loc[row, col] == data.loc[row, col]:
                failing.append('%s\t%s' % (samp, col))
            rginfo[samp][col] = data.loc[row, col]
        # optional RG info
        for col in ['rgid', 'rgpu']:
            if data.loc[row, col] != data.loc[row, col]:
                # if nan
                rginfo[samp][col] = None
                if samp not in warning:
                    warning.append(samp)
            else:
                rginfo[samp][col] = data.loc[row, col]

        # map between file and pool/samp
        for f in [data.loc[row, 'file_name_r1'], data.loc[row, 'file_name_r2']]:
            f2pool[f] = pool
            f2samp[op.join(pooldir, f)] = samp

    # handle --rm_paralogs, --translate, --rm_repeats
    for pool in uni(data['pool_name']):
        # handle translating stitched genome to unstitched positions
        pool2translate[pool] = handle_translate(translate, pool2translate, poolref[pool], data, pool)

        # handle removing SNPs from repeat regions
        pool2repeatsfile[pool] = handle_repeats(repeats, pool2repeatsfile, poolref[pool], data, pool)

        # handle removing paralogs
        pool2paralogfile[pool] = handle_paralogs(paralogs, pool2paralogfile, data, pool, parentdir)

    # handle fails for rm_repeats/translate/rm_paralogs
    handle_dict_fails(pool2repeatsfile, pool2translate, pool2paralogfile, repeats, translate, paralogs, data, parentdir)

    # RG info failing/warnings
    handle_rg_fails(failing, warning, parentdir, data)

    pkldump(pool2repeatsfile, op.join(parentdir, 'repeat_regions.pkl'))
    pkldump(pool2paralogfile, op.join(parentdir, 'paralog_snps.pkl'))
    pkldump(pool2translate, op.join(parentdir, 'translate_snps.pkl'))
    pkldump(rginfo, op.join(parentdir, 'rginfo.pkl'))
    pkldump(ploidy, op.join(parentdir, 'ploidy.pkl'))
    pkldump(f2samp, op.join(parentdir, 'f2samp.pkl'))
    pkldump(poolsamps, op.join(parentdir, 'poolsamps.pkl'))
    pkldump(poolref, op.join(parentdir, 'poolref.pkl'))
    pkldump(adaptors, op.join(parentdir, 'adaptors.pkl'))
    pkldump(samp2pool, op.join(parentdir, 'samp2pool.pkl'))
    return f2pool, poolref


def check_reqs(parentdir):
    """Check for assumed exports."""
    print(Bcolors.BOLD + '\nChecking for exported variables' + Bcolors.ENDC)
    variables = ['SLURM_ACCOUNT', 'SBATCH_ACCOUNT', 'SALLOC_ACCOUNT',
                 'VARSCAN_DIR', 'PYTHONPATH', 'SQUEUE_FORMAT']

    # check to see if bash_variables file has been created
    if not op.exists(op.join(parentdir, 'bash_variables')):
        print('\tCould not find bash_variables file in parentdir. Please create this file and add \
in variables from README (eg SLURM_ACCOUNT, SQUEUE_FORMAT, etc). See example in $HOME/pipeline.')
        print('exiting pipeline')
        exit()
    else:
        with open(op.join(parentdir, 'bash_variables')) as bv:
            text = bv.read().split("\n")
        needed = []
        for var in variables:
            found = False
            for line in text:
                if var in line:
                    found = True
                    break
            if found is False:
                needed.append(var)
        if len(needed) > 0:
            print(Bcolors.FAIL + '\tFAIL: not all bash variables were found in parentdir/bash_variables file.' + Bcolors.ENDC)
            print(Bcolors.FAIL + '\tFAIL: the following variables must be present' + Bcolors.ENDC)
            for var in needed:
                print(Bcolors.FAIL + '\t%s' % var + Bcolors.ENDC)
            print('exiting pipeline')

    # check to see if bash_variables file has been sourced
    for var in variables:
        try:
            print('\t%s = %s' % (var, os.environ[var]))
        except KeyError:
            print(Bcolors.FAIL + '\tCould not find %s in exported vars\n\texport this var in parentdir/bash_variables \
so it can be used later in pipeline, then source this file before restarting pipeline.' % var + Bcolors.ENDC)
            print('\texiting 00_start-pipeline.py')
            exit()

    # check for programs
    for program in [op.join(os.environ['VARSCAN_DIR'], 'VarScan.v2.4.3.jar')]:
        if not op.exists(program):
            print(Bcolors.BOLD +
                  Bcolors.FAIL +
                  "FAIL: could not find the following program: %s" % program +
                  Bcolors.ENDC)

    # make sure an environment can be activated (activation assumed to be in $HOME/.bashrc)
    for exe in ['activate']:
        if distutils.spawn.find_executable(exe) is None:
            print('\tcould not find %s in $PATH\nexiting 00_start-pipeline.py' % exe)
            if exe == 'activate':
                print('\t\t(the lack of activate means that the python env is not correctly installed)')
            exit()
    # make sure pipeline can be accessed via $HOME/pipeline
    if not op.exists(op.join(os.environ['HOME'], 'pipeline')):
        print('\tcould not find pipeline via $HOME/pipeline\n\texiting 00_start-pipeline.py')
        exit()


def check_pyversion():
    """Make sure python version is 3.6+"""
    pyversion = float(str(sys.version_info[0]) + '.' + str(sys.version_info[1]))
    if not pyversion >= 3.6:
        text = '''FAIL: You are using python %s. This pipeline was built with python 3.7+.
FAIL: You will need at least python v3.6+.
FAIL: exiting 00_start-pipeline.py
    ''' % pyversion
        print(Bcolors.BOLD + Bcolors.FAIL + text + Bcolors.ENDC)
        exit()


def get_pars():
    choices = ['all', 'fail', 'begin', 'end', 'pipeline-finish']
    parser = argparse.ArgumentParser(description=mytext,
                                     add_help=False,
                                     formatter_class=argparse.RawTextHelpFormatter)
    requiredNAMED = parser.add_argument_group('required arguments')
    requiredNAMED.add_argument("-p",
                               required=True,
                               default=argparse.SUPPRESS,
                               dest="parentdir",
                               type=str,
                               help="/path/to/directory/with/fastq.gz-files/")
    parser.add_argument("-e",
                        required=False,
                        dest="email",
                        help='''the email address you would like to have notifications
sent to''')
    parser.add_argument("-n",
                        default=None,
                        nargs='+',
                        required=False,
                        dest="email_options",
                        help='''the type(s) of email notifications you would like to
receive from the pipeline. Requires --email-address.
These options are used to fill out the #SBATCH flags.
Must be one (or multiple) of 
%s
(default: None)''' % [x for x in choices])
    parser.add_argument("-maf",
                        required=False,
                        dest="maf",
                        help='''At the end of the pipeline, VCF files will be filtered
for MAF. If the pipeline is run on a single
population/pool, the user can set MAF to 0.0 so as to
filter variants based on global allele frequency 
across populations/pools at a later time. (if the
number of sample_names in a pool == 1 then default
maf=0; Otherwise default maf = 1/sum(ploidy column)''')
    parser.add_argument('--translate',
                        required=False,
                        action='store_true',
                        dest="translate",
                        help='''Boolean: true if used, false otherwise. If a stitched
genome is used for mapping, this option will look for
a ref.order file in the same directory as the
ref.fasta - where ref is the basename of the ref.fasta
(without the .fasta). The pipeline will use this
.order file to translate mapped positions to
unstitched positions at the end of the pipeline while
filtering. Positions in .order file are assumed to be
1-based indexing. Assumes .order file has no header,
and is of the format (contig name from unstitched
genome, start/stop are positions in the stitched genome):
ref_scaffold<tab>contig_name<tab>start_pos<tab>stop_pos<tab>contig_length
(default: False)''')
    parser.add_argument('--rm_repeats',
                        required=False,
                        action='store_true',
                        dest='repeats',
                        help='''Boolean: true if used, false otherwise. If repeat
regions are available, remove SNPs that fall within
these regions from final SNP table and write to 
a REPEATS table. This option will look for a .txt file
in the same directory as the ref.fasta. Assumes the
filename is of the form: ref_repeats.txt - where ref
is the basename of the ref.fasta (without the .fasta).
This file should have 1-based indexing and should be
located in the same directory as the reference. The
file should have a header ('CHROM', 'start', 'stop').
The CHROM column can be names in the reference (if
using unstitched reference), or names of contigs that
were stitched to form the reference. If using a
stitched genome, --translate is required. (default:
False)''')
    parser.add_argument('--rm_paralogs',
                        required=False,
                        action='store_true',
                        dest='paralogs',
                        help='''Boolean: true if used, false otherwise. If candidate
sites have been isolated within the reference where
distinct gene copies (paralogs) map to the same
position (and thus create erroneous SNPs), remove any
SNPs that fall on these exact sites and write to a
PARALOGS file. The pipeline assumes this file is
located in the parentdir, andends with 
'_paralog_snps.txt'. This file is tab-delimited, and
must have a column called 'locus' thatcontains
hyphen-separated CHROM-POS sites for paralogs. These
sites should be found in the current ref.fa being
used to call SNPs (otherwise SNPs cannot be filtered
by these sites). (default: False)''')
    parser.add_argument('-h', '--help',
                        action='help',
                        default=argparse.SUPPRESS,
                        help='Show this help message and exit.\n')
    args = parser.parse_args()
    # trim path
    if args.parentdir.endswith('/'):
        args.parentdir = args.parentdir[:-1]
    # save command
    pkldump(args, op.join(args.parentdir, 'pipeline_start_command.pkl'))
    # assess arguments
    if args.email and args.email_options is None:
        print(Bcolors.FAIL + 'FAIL: --notification-types are required when specifying email' + Bcolors.ENDC)
        print(Bcolors.FAIL + 'FAIL: choices = {%s}\n' % [x for x in choices] + Bcolors.ENDC)
        exit()
    if args.email_options and args.email is None:
        print(Bcolors.FAIL + 'FAIL: specifying --notification-types requires specifying \
--email-address\n' + Bcolors.ENDC)
        exit()
    if args.email_options:
        for choice in args.email_options:
            if not choice.lower() in choices:
                print(Bcolors.FAIL +
                      '''FAIL: There can be multiple options, but they must be from the set:''' +
                      Bcolors.ENDC)
                print(Bcolors.FAIL +
                      '''\t%s\n''' % choices +
                      Bcolors.ENDC)
                exit()
    if args.email:
        if '@' not in args.email:
            print(Bcolors.FAIL + 'FAIL: email address does not have an "@" symbol in it, \
please check input\n' + Bcolors.ENDC)
            exit()
        if 'all' in args.email_options:
            args.email_options = ['all']
        # save email
        epkl = {'email': args.email,
                'opts': args.email_options}
        pkldump(epkl, op.join(args.parentdir, 'email_opts.pkl'))

    if args.maf:
        pkldump(args.maf, op.join(args.parentdir, 'maf.pkl'))

    if args.repeats:
        text = 'WARN: You have indicated that you want to remove repeats.\n'
        text = text + 'WARN: Make sure --translate is used if using a stitched reference.\n'
        text = text + 'WARN: Otherwise this will cause an error.\n'
        text = text + 'WARN: --repeats assumes that the first column in the repeats file ...\n'
        text = text + 'WARN: ... are the exact chromosome names found in the ref.fasta, ...\n'
        text = text + 'WARN: ... or if used with --translate this assumes that the first ...\n'
        text = text + 'WARN: ... column of the repeats file are names found in the second ...\n'
        text = text + 'WARN: ... column of the ref.order file used to translate positions.'
        print(Bcolors.WARNING + text + Bcolors.ENDC)
        askforinput()

    return args


def main():
    # parse arguments
    args = get_pars()

    # WARN if version = 3.6, FAIL if < 3.6
    check_pyversion()

    # look for exported vars (should be in .bashrc)
    check_reqs(args.parentdir)

    # determine which slurm accounts to use
    balance_queue.get_avail_accounts(args.parentdir, save=True)

    # read in the datatable
    data = read_datatable(args.parentdir)
    
    # create directories for each group of pools to be combined
    pooldirs = make_pooldirs(data, args.parentdir)
    
    # parse the datatable
    f2pool, poolref = parse_datatable(data,
                                      args.parentdir,
                                      args.translate,
                                      args.repeats,
                                      args.paralogs)

    # create bedfiles to parallelize varscan later on
    create_all_bedfiles(poolref, len(pooldirs))

    # create bedfiles to parallelize varscan later on
    create_all_bedfiles(poolref, len(pooldirs))

    # assign fq files to pooldirs for visualization (good to double check)
    get_datafiles(args.parentdir, f2pool, data)

    # create and sbatch sh files
    create_sh(pooldirs, poolref, args.parentdir)


if __name__ == '__main__':
    mytext = Bcolors.BOLD + Bcolors.OKGREEN + '''
*****************************************************************************


         ___|               \         |          _   __|
        |      _ \           \    __  |  _     _    |    _|  _ \\  _ \\
        |     (   | __|   /_  \  (    | (   | (  |  |   |    __/  __/
         ___|\___/      _/    _\\\___/_|\__/_|  __/  |  _|  \___|\___|
                                              |
                                              |

                             VarScan pipeline

*****************************************************************************


''' + Bcolors.ENDC
    main()
