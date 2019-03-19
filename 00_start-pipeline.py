"""
### fix
# uncomment create_sh.subprocess.call
# add in account argparse
###

### usage
# python 00_start-pipeline.py /path/to/folder/with/all/fastq/files/
###
"""

import os, sys, distutils.spawn, balance_queue, subprocess, shutil, argparse, create_bedfiles, pandas as pd
from os import path as op
from collections import OrderedDict
from coadaptree import fs, pkldump, uni, luni, makedir


def create_sh(pooldirs, poolref):
    # create sh files
    print('\nwriting sh files')
    for pooldir in pooldirs:
        pool = op.basename(pooldir)
        print('\npool = %s' % pool)
        ref = poolref[pool]
        print('sending pooldir and ref to 01_trim-fastq.py')
        #subprocess.call([shutil.which('python'), op.join(os.environ['HOME'], 'pipeline/01_trim-fastq.py'), pooldir, ref])
    print('\n')


def get_datafiles(parentdir, f2pool, data):
    # get list of files from datatable, make sure they exist in parentdir, create symlinks in /parentdir/<pool_name>/
    print('\nchecking for existance of fastq files in datatable.txt')
    files = [f for f in fs(parentdir) if 'fastq' in f and 'md5' not in f]
    datafiles = data['file_name_r1'].tolist()
    for x in data['file_name_r2'].tolist():
        datafiles.append(x)

    for f in datafiles:
        src = op.join(parentdir, f)
        if not op.exists(src)   :
            # make sure file in datatable exists
            print("could not find %s in %s\nmake sure file_name in datatable is its basename" % (f, parentdir))
            print("(symlinks in parentdir to fastq files in other dirs works fine, and is the intentional use)")
            sys.exit(1)
        pooldir = op.join(parentdir, f2pool[f])
        dst = op.join(pooldir, f)
        if not op.exists(dst):
            # easy to visualize in cmdline if script is finding correct group of files by ls-ing pooldir
            os.symlink(src, dst)


def make_pooldirs(data, parentdir):
    # make pool dirs
    print("\nmaking pool dirs")
    pools = uni(data['pool_name'].tolist())
    pooldirs = []
    for p in pools:
        DIR = op.join(parentdir, p)
        pooldirs.append(makedir(DIR))
    return pooldirs


def create_crisp_bedfiles(poolref):
    # create bedfiles for crisp
    print("\ncreating CRISP bedfiles")
    for ref in poolref.values():
        create_bedfiles.main('create_bedfiles.py', ref)


def read_datatable(parentdir):
    # read in the datatable, save rginfo for later
    datatable = op.join(parentdir, 'datatable.txt')
    if not op.exists(datatable):
        print('the datatable is not in the necessary path: %s\nexiting 00_start-pipeline.py' % datatable)
        sys.exit(3)
    print('reading datatable, getting fastq info')
    data = pd.read_csv(datatable, sep='\t')
    rginfo = {}     # key=sampname vals=rginfo
    samp2pool = {}  # key=samp val=pool
    poolref = {}    # key=pool val=ref.fa
    ploidy = {}     # key=pool val=ploidy
    poolsamps = {}  # key=pool val=sampnames
    f2samp = {}     # key=f val=samp
    f2pool = {}     # key=f val=pool
    adaptors = OrderedDict()  # key=samp val={'r1','r2'} val=adaptor
    for row in data.index:
        samp = data.loc[row, 'sample_name']
        adaptors[samp] = {'r1': data.loc[row, 'adaptor_1'],
                          'r2': data.loc[row, 'adaptor_2']}
        pool = data.loc[row, 'pool_name']
        pooldir = op.join(parentdir, pool)
        print('{}\tsamp = {}\tpool = {}'.format(row, samp, pool))
        if pool not in poolsamps:
            poolsamps[pool] = []
        if samp not in poolsamps[pool]:
            poolsamps[pool].append(samp)
        samp2pool[samp] = pool
        df = data[data['pool_name'] == pool].copy()
        if not luni(df['ploidy']) == 1:
            print("the ploidy values for some elements with pool name '%s' are not the same" % pool)
            sys.exit(1)
        if pool not in ploidy:
            ploidy[pool] = data.loc[row, 'ploidy']
        if pool in poolref:
            if not poolref[pool] == data.loc[row, 'ref']:
                print("ref genome for samples in %s pool seems to have different paths in datatable.txt" % pool)
                sys.exit(1)
        else:
            ref = data.loc[row, 'ref']
            if not op.exists(ref):
                print('ref for %s does not exist in path: %s' % (samp, ref))
                print('exiting %s' % thisfile)
                exit()
            poolref[pool] = ref
        rginfo[samp] = {}
        for col in ['rglb', 'rgpl', 'rgsm']:  # rg info columns
            rginfo[samp][col] = data.loc[row, col]
        for f in [data.loc[row, 'file_name_r1'], data.loc[row, 'file_name_r2']]:
            f2pool[f] = pool
            f2samp[op.join(pooldir, f)] = samp
    pkldump(rginfo, op.join(parentdir, 'rginfo.pkl'))
    pkldump(ploidy, op.join(parentdir, 'ploidy.pkl'))
    pkldump(f2samp, op.join(parentdir, 'f2samp.pkl'))
    pkldump(poolsamps, op.join(parentdir, 'poolsamps.pkl'))
    pkldump(poolref, op.join(parentdir, 'poolref.pkl'))
    pkldump(adaptors, op.join(parentdir, 'adaptors.pkl'))
    return data, f2pool, poolref


def check_reqs():
    # check for assumed exports
    print('\nchecking for exported variables')
    for var in ['SLURM_ACCOUNT', 'SBATCH_ACCOUNT', 'SALLOC_ACCOUNT',
                'CRISP_DIR', 'PYTHONPATH', 'SQUEUE_FORMAT']:
        try:
            print('\t%s = %s' % (var, os.environ[var]))
        except KeyError:
            print('\tcould not find %s in exported vars\n\texport this var in $HOME/.bashrc so it can be used later in pipeline\n\texiting 00_start-pipeline.py' % var)
            exit()
    for exe in ['lofreq', 'activate']:
        if distutils.spawn.find_executable(exe) is None:
            print('\tcould not find %s in $PATH\nexiting 00_start-pipeline.py' % exe)
            if exe == 'activate':
                print('\t\t(the lack of activate means that the python env is not correctly installed)')
            exit()
    print('DONE!\n')


def check_pyversion():
    # check python version
    pyversion = str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + '.' + str(sys.version_info[2])
    if not sys.version_info[0] == 3:
        text = '''FAIL: You are using python %s. This pipeline was built with python 3.7+.
    FAIL: Please use a more recent version of python.
    FAIL: exiting 00_start-pipeline.py
    ''' % pyversion
        print(Bcolors.BOLD + Bcolors.FAIL + text + Bcolors.ENDC)
        exit()
    if not sys.version_info[1] == 7:
        text = "WARN: You are using python v%s. This pipeline was built with python v3.7+.\n" \
               "WARN: You may want to consider updating to a more recent version of python.\n" % pyversion
        print(Bcolors.BOLD + Bcolors.WARNING + text + Bcolors.ENDC)
        while True:
            inp = input("INPUT NEEDED: Do you want to proceed? (yes | no): ").lower()
            if inp in ['yes', 'no']:
                break
            else:
                print("Please respond with 'yes' or 'no'")
        if inp == 'no':
            print('exiting 00_start-pipeline.py')
            exit()


def get_pars():
    parser = argparse.ArgumentParser(description=print(mytext),
                                     add_help=False,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-p", "--parent-dir",
                        required=True,
                        dest="parentdir",
                        type=str,
                        help="/path/to/directory/with/fastq.gz-and-datatabe.txt/files/")
    parser.add_argument("-e", "--email-address",
                        required=False,
                        dest="email",
                        help="the email address you would like to have notifications sent to")
    parser.add_argument("-n", "--notification-types",
                        default=['pipeline-finish'],
                        nargs='+',
                        choices=['all','none','fail','begin','end','pipeline-finish'],
                        required=False,
                        dest="email_options",
                        help="the type(s) of email notifications you would like to receive from the pipeline.\
                        Requires --email-address.")
    parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                        help='Show this help message and exit.\n')
    args = parser.parse_args()
    if args.parentdir.endswith('/'):
        args.parentdir = args.parentdir[:-1]
    if args.email and args.email_options is None:
        parser.error('--email-address requires --notification-types')
    if args.email_options and args.email is None:
        parser.error('specifying --notification-types requires specifying --email-address')
    if args.email:
        if not '@' in args.email:
            parser.error('email address does not have an "@" symbol in it, please check input')
        if 'all' in args.email_options:
            args.email_options = ['all']
        # save email
        epkl = {'email': args.email,
                'opts': args.email_options}
        pkldump(epkl,op.join(args.parentdir,'email_opts.pkl'))
    print(vars(args))
    print(args.parentdir)

    return args


def main():
    # parse arguments
    args = get_pars()

    # make sure version >= 3, WARN if < 3.7
    check_pyversion()

    # look for exported vars (should be in .bashrc)
    check_reqs()

    data, f2pool, poolref = read_datatable(args.parentdir)

    create_crisp_bedfiles(poolref)

    pooldirs = make_pooldirs(data, args.parentdir)

    get_datafiles(args.parentdir, f2pool, data)

    create_sh(pooldirs, poolref)


if __name__ == '__main__':
    class Bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'
    mytext = Bcolors.BOLD + Bcolors.OKGREEN + '''
*****************************************************************************


   ___|                 \          |              _   __|
  |       _ \            \     __  |   _      _      |      _|    _ \\    _ \\
  |      (   | __|    /_  \   (    |  (   |  (  |    |     |      __/    __/
   ___| \___/       _/    _\ \___/_| \__/_|   __/    |    _|    \___|  \___|
                                             |
                                             |

                          LoFreq and CRISP pipeline

*****************************************************************************


''' + Bcolors.ENDC

    main()