"""
start the pipeline.

### usage
# 00_start-pipeline.py -p PARENTDIR [-e EMAIL [-n EMAIL_OPTIONS]]
###

### assumes
# that samples duplicated in the 'sample_name' column have the same rglb, rgsm, and rgpl read groups
# that ploidy is the same (ie uses the exact same pop members) across resquencing of ...
#     individual pools (all resequencing must have same ploidy for any given sample)
###
"""

import os, sys, distutils.spawn, subprocess, shutil, argparse, pandas as pd
import balance_queue, create_bedfiles
from os import path as op
from subprocess import PIPE
from subprocess import Popen
from collections import OrderedDict
from coadaptree import fs, pkldump, uni, makedir, askforinput, Bcolors

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
        askforinput()

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
        DIR = op.join(parentdir, p)
        if op.exists(DIR):
            print("The pooldir already exists, this could overwrite previous data: %s" % DIR)
            askforinput()
        pooldirs.append(makedir(DIR))
    return pooldirs


def create_all_bedfiles(poolref):
    """For each unique ref.fa in datatable.txt, create bedfiles for varscan.

    Positional arguments:
    poolref - dictionary with key = pool, val = /path/to/ref
    """
    # create bedfiles for varscan
    print(Bcolors.BOLD + "\ncreating bedfiles" + Bcolors.ENDC)
    for ref in uni(poolref.values()):
        create_bedfiles.main(ref)


def read_datatable(parentdir, translate, repeats, paralogs):
    """
    Read in datatable.txt, use to create files and dirs for downstream.
    Also checks some assumptions of datatable.txt.

    translate, repeats, and paralogs are boolean.
    parentdir is a path.
    """
    datatable = op.join(parentdir, 'datatable.txt')
    if not op.exists(datatable):
        print(Bcolors.FAIL + '''FAIL: the datatable is not in the necessary path: %s
FAIL: exiting 00_start-pipeline.py''' % datatable + Bcolors.ENDC)
        sys.exit(3)
    print(Bcolors.BOLD + 'reading datatable, getting fastq info' + Bcolors.ENDC)
    data = pd.read_csv(datatable, sep='\t')
    rginfo = {}     # key=samp vals=rginfo
    samp2pool = {}  # key=samp val=pool
    poolref = {}    # key=pool val=ref.fa
    ploidy = {}     # key=pool val=dict(key=sample: val=sample_ploidy)
    poolsamps = {}  # key=pool val=sampnames
    f2samp = {}     # key=f val=samp
    f2pool = {}     # key=f val=pool
    adaptors = OrderedDict()  # key=samp val={'r1','r2'} val=adaptor
    warning = [] # whether to print out warning about RGID
    for row in data.index:
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
        if pool not in ploidy:
            ploidy[pool] = {}
        if samp in ploidy[pool].keys():
            if ploidy[pool][samp] != int(data.loc[row, 'ploidy']):
                text = "FAIL: the ploidy values for sample_name '%s' are not the same" % samp
                print(Bcolors.FAIL + text + Bcolors.ENDC)
                exit()
        ploidy[pool][samp] = int(data.loc[row, 'ploidy'])
        if pool in poolref:
            if not poolref[pool] == data.loc[row, 'ref']:
                print("ref genome for samples in %s pool seems to have different paths in datatable.txt" % pool)
                sys.exit(1)
        else:
            ref = data.loc[row, 'ref']
            if not op.exists(ref):
                text = 'FAIL: ref for %s does not exist in path: %s' % (samp, ref)
                print(Bcolors.FAIL + text + Bcolors.ENDC)
                print('exiting 00_start-pipeline.py')
                exit()
            if translate is True:
                orderfile = ref.split(".fa")[0] + '.order'
                if not op.exists(orderfile):
                    text = 'FAIL: You have indicated that you would like stitched regions translated.'
                    text = text + 'But the pipeline cannot find the .order file: %s' % orderfile
                    text = text + '\nexiting 00_start-pipeline.py'
                    print(Bcolors.FAIL + text + Bcolors.ENDC)
                    exit()
                else:
                    pkldump(orderfile, op.join(parentdir, 'orderfile.pkl'))
            if repeats is True:
                repeatfile = ref.split(".fa")[0] + '_repeats.txt'
                if not op.exists(repeatfile):
                    text = 'FAIL: You have indicated that you would like repeat regions removed. '
                    text = text + 'But the pipeline cannot find the file for repeat regions: %s' % repeatfile
                    print(Bcolors.FAIL + text + Bcolors.ENDC)
                    exit()
                else:
                    pkldump(repeatfile, op.join(parentdir, 'repeat_regions.pkl'))
            if paralogs is True:
                parafiles = [f for f in fs(parentdir) if f.endswith('_paralog_snps.txt')]
                if len(parafiles) > 1:
                    text = f'FAIL: There are multiple files in {parentdir} with "_paralog_snps.txt" in the name.\n'
                    text = text + 'FAIL: Please remove all but one of these files.'
                    print(Bcolors.FAIL + text + Bcolors.ENDC)
                    print('exiting 00_start.py')
                    exit()
                elif len(parafiles) == 0:
                    parafile = ''
                elif len(parafiles) == 1:
                    parafile = parafiles[0]
                if not op.exists(parafile):
                    text = 'FAIL: You have indicated that you would like paralog sites removed '
                    text = text + 'from final SNPs. But the pipeline cannot find the file for repeat '
                    text = text + 'regions. This file must end with "_paralog_snps.txt".'
                    print(Bcolors.FAIL + text + Bcolors.ENDC)
                    exit()
                else:
                    pkldump(parafile, op.join(parentdir, 'paralog_snps.pkl'))
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
            poolref[pool] = ref
        rginfo[samp] = {}
        for col in ['rglb', 'rgpl', 'rgsm']:  # rg info columns
            rginfo[samp][col] = data.loc[row, col]
        for col in ['rgid', 'rgpu']:
            if data.loc[row, col] != data.loc[row, col]:
                # if nan
                rginfo[samp][col] = None
                if samp not in warning:
                    warning.append(samp)
            else:
                rginfo[samp][col] = data.loc[row, col]
        for f in [data.loc[row, 'file_name_r1'], data.loc[row, 'file_name_r2']]:
            f2pool[f] = pool
            f2samp[op.join(pooldir, f)] = samp
    if len(warning) > 0:
        outputs = []
        for row in data.index:
            samp = data.loc[row, 'sample_name']
            if samp in warning:
                r1 = op.join(parentdir, data.loc[row, 'file_name_r1'])
                outputs.append("%s\t%s" % (samp, get_rgid(r1)))
        print(Bcolors.WARNING + 'WARN: at least one of the samples has a blank RGID.\n' +
              'WARN: The pipeline will automatically assign the following RGIDs.\n' +
              Bcolors.ENDC)
        for output in outputs:
            print(Bcolors.WARNING + output + Bcolors.ENDC)
        print('\n', Bcolors.WARNING + 'WARN: If RGPU is also blank, the pipeline will assign RGPU as: $RGID.$RGLB' +
              Bcolors.ENDC)
        askforinput()
                            
                               
    pkldump(rginfo, op.join(parentdir, 'rginfo.pkl'))
    pkldump(ploidy, op.join(parentdir, 'ploidy.pkl'))
    pkldump(f2samp, op.join(parentdir, 'f2samp.pkl'))
    pkldump(poolsamps, op.join(parentdir, 'poolsamps.pkl'))
    pkldump(poolref, op.join(parentdir, 'poolref.pkl'))
    pkldump(adaptors, op.join(parentdir, 'adaptors.pkl'))
    pkldump(samp2pool, op.join(parentdir, 'samp2pool.pkl'))
    return data, f2pool, poolref


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
    parser = argparse.ArgumentParser(description=print(mytext),
                                     add_help=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
                        help="the email address you would like to have notifications sent to")
    parser.add_argument("-n",
                        default=None,
                        nargs='+',
                        required=False,
                        dest="email_options",
                        help='''the type(s) of email notifications you would like to receive from the pipeline.\
                        Requires --email-address. These options are used to fill out the #SBATCH flags.
must be one (or multiple) of %s''' % [x for x in choices])
    parser.add_argument("-maf",
                        required=False,
                        dest="maf",
                        help='''At the end of the pipeline, VCF files will be filtered for MAF. If the pipeline is run on a single population/pool, the user can set MAF to 0.0 so as to filter variants based on global allele frequency across populations/pools at a later time. (if the number of sample_names in a pool == 1 then default maf=0; Otherwise maf = 1/sum(ploidy column)''')
    parser.add_argument('--translate',
                        required=False,
                        action='store_true',
                        dest="translate",
                        help='''Boolean: true if used, false otherwise. If a stitched genome is used for mapping, this option will look for a ref.order file in the same directory as the ref.fasta - where ref is the basename of the ref.fasta (without the .fasta). The pipeline will use this .order file to translate mapped positions to unstitched positions at the end of the pipeline while filtering. Positions in .order file are assumed to be 1-based indexing. Assumes .order file has no header, and is of the format (contig name from unstitched genome, start/stop are positions in the stitched genome):
ref_scaffold<tab>contig_name<tab>start_pos<tab>stop_pos<tab>contig_length''')
    parser.add_argument('--rm_repeats',
                        required=False,
                        action='store_true',
                        dest='repeats',
                        help="Boolean: true if used, false otherwise. If repeat regions are available, remove SNPs that fall within these regions. This option will look for a .txt file in the same directory as the ref.fasta. Assumes the filename is of the form: ref_repeats.txt - where ref is the basename of the ref.fasta (without the .fasta). This file should have 1-based indexing and should be located in the same directory as the reference. The file should have a header ('CHROM', 'start', 'stop'). The CHROM column can be names in the reference (if using unstitched reference), or names of contigs that were stitched to form the reference. If using a stitched genome, --translate is required.")
    parser.add_argument('--rm_paralogs',
                        required=False,
                        action='store_true',
                        dest='paralogs',
                        help="Boolean: true if used, false otherwise. If candidate sites have been isolated within the reference where distinct gene copies (paralogs) map to the same position (and thus create erroneous SNPs), remove any SNPs that fall on these exact sites. The pipeline assumes this file is located in the parentdir, and ends with '_paralog_snps.txt'. This file is tab-delimited, and must have a column called 'locus' that contains hyphen-separated CHROM-POS sites for paralogs. These sites should be found in the current ref.fa being used to call SNPs (otherwise SNPs cannot be filtered by these sites).")
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
    data, f2pool, poolref = read_datatable(args.parentdir,
                                           args.translate,
                                           args.repeats,
                                           args.paralogs)

    # create bedfiles to parallelize varscan later on
    create_all_bedfiles(poolref)

    # create directories for each group of pools to be combined
    pooldirs = make_pooldirs(data, args.parentdir)

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
