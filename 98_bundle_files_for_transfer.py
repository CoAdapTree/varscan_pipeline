"""Create scp commands to copy completed files from current dir to remote server\
(from perspective of remote server).

# usage
# python 98_bundle_files_for_transfer.py parentdir remote_parentdir <bool>
# 

# assumes
# that transfer is done in parallel on remote server
# that compute canada servers are abbreviated in your remote:$HOME/.ssh/config (eg cedar, beluga, graham)
# that any md5 files are for the current files in the directory that also haven't been modified since md5 creation
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
    """See if an .md5 file exists, if not create .md5 file.
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
    subcmds = []
    for src in srcfiles:
        if createmd5 is True:
            md5 = check_md5(src, md5files)
            md5dst = op.join(remotedir, op.basename(md5))
            subcmds.append(f'scp {hostname}:{md5} {md5dst}')
        dst = op.join(remotedir, op.basename(src))
        subcmds.append(f'scp {hostname}:{src} {dst}')
    return subcmds


pools = list(pklload(op.join(parentdir, 'poolref.pkl')).keys())
pooldirs = [op.join(parentdir, p) for p in pools]
newdirs = []  # keep track of directories to easily make on remote server
cmds = []  # keep track of all scp commands
# get hostname (eg from beluga2.int.ets1.calculquebec.ca)
hostname = ''.join([x for x in os.environ['HOSTNAME'].split(".")[0] if x.isalpha()])


# get shfiles
for p in pooldirs:
    shdir = op.join(p, 'shfiles')
    remotesh = op.join(remote, f'{p}/sh_and_outfiles')
    newdirs.append(remotesh)
    dirs = [d for d in fs(shdir) if op.isdir(d)]
    for d in dirs:
        remoted = op.join(remotesh, op.basename(d))
        newdirs.append(remoted)
        md5files = [f for f in fs(d) if f.endswith('.md5')]
        srcfiles = [f for f in fs(d) if not f.endswith('.md5')]
        cmds.extend(get_cmds(srcfiles, md5files, remoted, False))


# get realigned bamfiles
for p in pooldirs:
    bamdir = op.join(p, '04_realign')
    remotebamdir = op.join(remote, f'{p}/realigned_bamfiles')
    newdirs.append(remotebamdir)
    md5files = [f for f in fs(bamdir) if f.endswith('.md5')]
    srcfiles = [f for f in fs(bamdir) if not f.endswith('.md5')]
    cmds.extend(get_cmds(srcfiles, md5files, remotebamdir, generate_md5))


# get read info
readinfo = op.join(parentdir, 'readinfo.txt')
readinfodst = op.join(remote, 'readinfo.txt')  # no need to add to newdirs
if not op.exists(readinfo):
    warning = "WARN: readinfo.txt does not exist. (you can run 99_get_read_stats.py and transfer later)"
    print(Bcolors.WARNING + warning + Bcolors.ENDC)
    askforinput()
else:
    cmds.append(f"scp {hostname}:{readinfo} {readinfodst}")


# get datatable used
datatable = op.join(parentdir, 'datatable.txt')
datatabledst = op.join(remote, 'datatable.txt')
cmds.append(f"scp {hostname}:{datatable} {datatabledst}")


# get varscan output
for p in pooldirs:
    varscan = op.join(p, 'varscan')
    remotevarscan = op.join(remote, f'{p}/snpsANDindels')
    newdirs.append(remotevarscan)
    md5files = [f for f in fs(varscan) if f.endswith('.md5') and '_all_' not in f]
    srcfiles = [f for f in fs(varscan) if f.endswith('.gz') or f.endswith('.txt') and '_all_' not in f]
    cmds.extend(get_cmds(srcfiles, md5files, remotevarscan, generate_md5))
    # double check for _all_SNPs and _all_INDELs
    md5files = [f for f in fs(varscan) if f.endswith('.md5') and '_all_' in f]
    srcfiles = [f for f in fs(varscan) if f.endswith('.txt') and '_all_' in f]
    if not len(srcfiles) == 2:
        warning = "\nWARN: There are not two all-files (SNP + INDEL) which are expected output"
        warning = warning + "\nWARN: Here are the files I found:\n"
        warning = warning + "\n\t".join(srcfiles)
        print(Bcolors.BOLD + Bcolors.WARNING + warning + Bcolors.ENDC)
        askforinput()
    cmds.extend(get_cmds(srcfiles, md5files, remotevarscan, generate_md5))


# write commands to file
scpfile = op.join(parentdir, 'scp_cmds.txt')
print(Bcolors.BOLD + f'\nwriting {len(cmds)} commands to: ' + Bcolors.ENDC + f'{scpfile}')
with open(scpfile, 'w') as o:
    jcmds = '\n'.join(cmds)
    o.write("%s" % jcmds)


# print out necessary dirs on remote
print(Bcolors.BOLD + "\nCreate the following dirs on remote before executing scp commands:" + Bcolors.ENDC)
for d in newdirs:
    print('mkdir ', d)
