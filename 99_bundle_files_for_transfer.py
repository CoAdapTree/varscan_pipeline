"""Create rsync commands to rsync completed files from current dir to remote server \
(from perspective of remote server).

# usage of this script
# python 98_bundle_files_for_transfer.py parentdir remote_parentdir <bool>
#

# usage of rsync_cmds.txt file output from this script
# [remote_server]$ cat rsync_cmds.txt | parallel -j 56 --progress

# assumes
# that transfer is done in parallel on remote server
# that compute canada servers are abbreviated in your remote:$HOME/.ssh/config
#    (eg cedar, beluga, graham)
# that any md5 files are for the current files in the directory that also haven't
#    been modified since md5 creation

"""

import os, sys, subprocess, shutil
from os import path as op
from coadaptree import fs, askforinput, Bcolors, pklload


thisfile, parentdir, remote, generate_md5 = sys.argv
generate_md5 = True if generate_md5 == 'True' else False
if parentdir.endswith("/"):
    parentdir = parendir[:-1]
if remote.endswith("/"):
    remote = remote[:-1]
print(Bcolors.BOLD + "input arguments:" + Bcolors.ENDC)
print('\t', 'parentdir = ', parentdir)
print('\t', 'remote = ', remote)
print('\t', 'generate_md5 = ', generate_md5)


def check_md5(src, md5files):
    """
    See if an .md5 file exists, if not create .md5 file.
    Only used for non-sh/out files.
    """
    os.chdir(op.dirname(src))
    md5 = src + '.md5'
    if md5 not in md5files and generate_md5 is True:
        print(f'creating md5 for {op.basename(src)} ...')
        md5hash = subprocess.check_output([shutil.which('md5sum'),
                                           src]).decode('utf-8').split()[0]
        with open(md5, 'w') as o:
            o.write(f"{md5hash}  {op.basename(src)}")
    return md5


def get_cmds(srcfiles, md5files, remotedir, createmd5):
    """
    Create rsync command between src_server and remote_server.
    Create .md5 if necessary.
    """
    subcmds = []
    for src in srcfiles:
        if createmd5 is True:
            md5 = check_md5(src, md5files)
            md5dst = op.join(remotedir, op.basename(md5))
            subcmds.append(f'rsync -azv {hostname}:{md5} {md5dst}')
        dst = op.join(remotedir, op.basename(src))
        subcmds.append(f'rsync -azv {hostname}:{src} {dst}')
    return subcmds


pools = list(pklload(op.join(parentdir, 'poolref.pkl')).keys())
pooldirs = [op.join(parentdir, p) for p in pools]
newdirs = []  # keep track of directories to easily make on remote server
cmds = []  # keep track of all rsync commands
# get hostname (eg beluga, cedar, graham)
hostname = os.environ['CC_CLUSTER']


# add remote and subdirs to newdirs list
newdirs.append(remote)
for p in pooldirs:
    newdirs.append(op.join(remote, op.basename(p)))


# get pkl files
print(Bcolors.BOLD + '\nBundling .pkl files ...' + Bcolors.ENDC)
pkls = [f for f in fs(parentdir) if f.endswith('.pkl')]
for p in pooldirs:
    for pkl in pkls:
        pkldst = op.join(remote, f'{op.basename(p)}/{op.basename(pkl)}')
        cmds.append(f"rsync -azv {hostname}:{pkl} {pkldst}")
    newpkls = [f for f in fs(p) if f.endswith('.pkl')]
    for newpkl in newpkls:
        newdst = op.join(remote, f'{op.basename(p)}/{op.basename(newpkl)}')
        cmds.append(f"rsync -azv {hostname}:{newpkl} {newdst}")


# get shfiles
print(Bcolors.BOLD + '\nBundling .sh and .out files ...' + Bcolors.ENDC)
for p in pooldirs:
    shdir = op.join(p, 'shfiles')
    remotesh = op.join(remote, f'{op.basename(p)}/sh_and_outfiles')
    newdirs.append(remotesh)
    dirs = [d for d in fs(shdir) if op.isdir(d)]
    for d in dirs:
        remoted = op.join(remotesh, op.basename(d))
        newdirs.append(remoted)
        md5files = [f for f in fs(d) if f.endswith('.md5')]
        srcfiles = [f for f in fs(d) if not f.endswith('.md5')]
        cmds.extend(get_cmds(srcfiles, md5files, remoted, False))


# get bedtools coords and samtools flagstat
print(Bcolors.BOLD + '\nBundling bedtools coords and samtools flagstats ...' + Bcolors.ENDC)
for p in pooldirs:
    bwadir = op.join(p, '02c_sorted_bamfiles')
    remotebwadir = op.join(remote, f'{op.basename(p)}/bedcoords_samflagstats')
    newdirs.append(remotebwadir)
    coords = [f for f in fs(bwadir) if 'coord' in f]
    flags = [f for f in fs(bwadir) if 'flagstat' in f]
    cmds.extend(get_cmds(coords, [], remotebwadir, False))
    cmds.extend(get_cmds(flags, [], remotebwadir, False))


# get realigned bamfiles
print(Bcolors.BOLD + '\nBundling realigned bamfiles ...' + Bcolors.ENDC)
for p in pooldirs:
    bamdir = op.join(p, '04_realign')
    remotebamdir = op.join(remote, f'{op.basename(p)}/realigned_bamfiles')
    newdirs.append(remotebamdir)
    md5files = [f for f in fs(bamdir) if f.endswith('.md5')]
    srcfiles = [f for f in fs(bamdir) if not f.endswith('.md5')]
    cmds.extend(get_cmds(srcfiles, md5files, remotebamdir, generate_md5))


# get read info
print(Bcolors.BOLD + '\nBundling readinfo.txt ...' + Bcolors.ENDC)
readinfo = op.join(parentdir, 'readinfo.txt')
datatable = op.join(parentdir, 'datatable.txt')
if not op.exists(readinfo):
    warning = "\tWARN: readinfo.txt does not exist. (you can run 99_get_read_stats.py and transfer later)"
    print(Bcolors.WARNING + warning + Bcolors.ENDC)
    askforinput(tab='\t', newline='')
else:
    for p in pooldirs:
        remotep = op.join(remote, op.basename(p))
        readinfodst = op.join(remotep, 'readinfo.txt')
        datatabledst = op.join(remotep, 'datatable.txt')
        cmds.append(f"rsync -azv {hostname}:{readinfo} {readinfodst}")  # no need to add to newdirs
        cmds.append(f"rsync -azv {hostname}:{datatable} {datatabledst}")  # no need to add to newdirs


# get varscan output
print(Bcolors.BOLD + '\nBundling varscan output ...' + Bcolors.ENDC)
paralogs = pklload(op.join(parentdir, 'paralog_snps.pkl'))  # used to determine if PARALOGS file is expected
repeats = pklload(op.join(parentdir, 'repeat_regions.pkl'))  # used to determine if REPEATS file is expected
poolseqcmd = vars(pklload(op.join(parentdir, 'pipeline_start_command.pkl')))
for p in pooldirs:
    pool = op.basename(p)
    varscan = op.join(p, 'varscan')
    if not op.exists(varscan):
        warning = f"\n\tWARN: varscan dir does not exist for pool: {op.basename(p)}"
        print(Bcolors.BOLD + Bcolors.WARNING + warning + Bcolors.ENDC)
        askforinput(tab='\t', newline='')
        continue
    remotevarscan = op.join(remote, f'{op.basename(p)}/snpsANDindels')
    newdirs.append(remotevarscan)
    remote_unfiltered = op.join(remotevarscan, '01_unfiltered')
    newdirs.append(remote_unfiltered)
    md5files = [f for f in fs(varscan) if f.endswith('.md5') and '_all_' not in f]
    srcfiles = [f for f in fs(varscan) if f.endswith('.gz') or f.endswith('.txt') and '_all_' not in f]
    cmds.extend(get_cmds(srcfiles, md5files, remote_unfiltered, generate_md5))
    # double check for _all_SNPs and _all_INDELs +/- _all_PARALOGS _all_REPEATS (baseline filtered)
    md5files = [f for f in fs(varscan) if f.endswith('.md5') and '_all_' in f]
    srcfiles = [f for f in fs(varscan) if f.endswith('.txt') and '_all_' in f]
    # determine the number of srcfiles that should be expected
    expected = ['SNP','INDEL']
    if poolseqcmd['repeats'] is True and repeats[pool] is not None:
        expected.append('REPEATS')
    if poolseqcmd['paralogs'] is True and paralogs[pool] is not None:
        expected.append('PARALOGS')
    if not len(srcfiles) == len(expected):
        warning = f"\n\tWARN: There are not {len(expected)} all-files ({' + '.join(expected)}) which are expected output for pool: {op.basename(p)}"
        warning = warning + f"\n\tWARN: Here are the {len(srcfiles)} files I found:\n"
        warning = warning + "\t" + "\n\t".join(srcfiles)
        print(Bcolors.BOLD + Bcolors.WARNING + warning + Bcolors.ENDC)
        askforinput(tab='\t', newline='')
    remote_filtered = op.join(remotevarscan, '02_baseline_filtered')
    newdirs.append(remote_filtered)
    cmds.extend(get_cmds(srcfiles, md5files, remote_filtered, generate_md5))


# write commands to file
print(Bcolors.BOLD + '\nWriting commands to rsync_cmds.txt file ...' + Bcolors.ENDC)
rsyncfile = op.join(parentdir, 'rsync_cmds.txt')
print(Bcolors.BOLD + f'\nwriting {len(cmds)} commands to: ' + Bcolors.ENDC + f'{rsyncfile}')
with open(rsyncfile, 'w') as o:
    jcmds = '\n'.join(cmds)
    o.write("%s" % jcmds)


# print out necessary dirs on remote
text = "\nCreate the following dirs on remote before executing rsync commands:"
text = text + "\n\t(or use find/replace in rsync file above)"
print(Bcolors.BOLD + text + Bcolors.ENDC)
for d in sorted(newdirs):
    print('mkdir ', d)
