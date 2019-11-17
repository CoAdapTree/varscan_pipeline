"""Create bedfiles for use in varscan and crisp.

### purpose
# create bedfiles from reference so that we can parallelize CRISP
# if an intervals directory exists, use .list files to create bedfiles (for stitched refs)
#    ELSE: look for .order file, and create bedfiles from this (for stitched refs)
#    ELSE: look for a .length file (create from ref if none exists) - assumes non-stitched
#

### usage
# python create_bedfiles.py /path/to/reference.fasta
###

### fix
# put bedfiles in parentdir instead of ref dir, so that all beds are the ones that are wanted
# make sure befile name str().zfill() matches digits of number used to define thresh
    # reduce duplicate text for calling thresh
###
"""

import sys, os, math
from os import path as op
from coadaptree import fs, makedir, askforinput, Bcolors


def openlenfile(lenfile):
    """
    Open lenfile to determine length of each contig in ref.fa.

    Positional arguments:
    lenfile = path to a file created from the ref.fa that has lengths of each contig/chromosome
    """
    with open(lenfile, 'r') as o:
        text = o.read().split("\n")
    return text


def make_beddir():
    """Create dir for bedfiles."""
    bname = op.basename(ref).split(".fa")[0]
    beddir = makedir(op.join(op.dirname(ref), 'bedfiles_%s' % bname))
    return bname, beddir


def get_prereqs(num):
    """Create a name for a bedfile based on the ref.fa path name and num.

    Positional arguments:
    num - int; the num'th bedfile
    """
    bname, beddir = make_beddir()
    f = op.join(beddir, "%s_bedfile_%s.bed" % (bname, str(num).zfill(4)))
    return f


def make_bed(lines, num):
    """Write contig/chrom, start, stop positions to bedfile.
    Different than make_bedfile(): .list files use zero-based, no need to correct."""
    f = get_prereqs(num)
    with open(f, 'w') as o:
        for contig, start, stop in lines:
            o.write("%s\t%s\t%s\n" % (contig, start, stop))


def make_bed_from_intervals(intdir):
    """If intervals.list files exist, use these instead of ref.fa.length file.

    Positional arguments:
    intdir - path to intervals.list files
    """
    intfiles = [f for f in fs(intdir) if f.endswith('.list')]
    for intfile in intfiles:
        num = intfile.split("_")[-1].replace(".list", "")
        lines = []
        with open(intfile, 'r') as o:
            text = o.read().split("\n")
        for line in text:
            scaff, span = line.split(":")
            start, stop = span.split("-")
            start, stop = (int(start) - 1, int(stop) - 1)
            lines.append((scaff, start, stop))
        make_bed(lines, num)
    print('\t\tcreated %s bedfiles for %s from interval files' % (len(intfiles), ref))


def make_beds_from_orderfile():
    """Use ref.order file to create bedfiles for parallelization."""
    orderfile = ref.replace(".fa", "") + '.order'
    print('\n\tCreating bedfiles from %s. Please confirm:\n\tAssuming .order file is of format:\n\t\tref_scaff<tab>contig_name<tab>start_pos<tab>stop_pos<tab>contig_length' % orderfile)
    askforinput()
    with open(orderfile, 'r') as o:
        text = o.read().split("\n")
    thresh = math.ceil(len(text) / 500)
    lines = []
    fcount = 0
    for count, line in enumerate(text):
        if not line == '':
            splits = line.split("\t")
            lines.append(([splits[0], int(splits[2]) - 1, int(splits[3]) - 1]))  # scaff start stop, zero-based
        if len(lines) == thresh or (count + 1 == len(text)):
            make_bedfile(lines, fcount, from_orderfile=True)
            lines = []
            fcount += 1


def find_positions():
    """Find positions to create bedfiles. First look for an intervals directory,
    then for a ref.order file. If neither exist, use ref.length file (ie, a lenfile)
    - if this doesn't exist, create one from the ref.fa.
    A lenfile is a file created from the ref.fa that has lengths of each contig/chromosome.
    """
    refdir = op.dirname(ref)
    intdir = op.join(refdir, 'intervals')
    # look for intervals directory with .list files
    if op.exists(intdir):
        print('\tusing intervals dir to create bedfiles for %s' % ref)
        make_bed_from_intervals(intdir)
        return
    # look for a ref.order file (assumed format will print)
    orderfile = ref.replace(".fa", "") + '.order'
    if op.exists(orderfile):
        make_beds_from_orderfile()
        return
    # create ref.fa.length if it doesn't exist
    if not op.exists('%(ref)s.length' % globals()):
        print("\tcreating %s.length file (this will take a few minutes)" % op.basename(ref))
        os.system('''cat  %(ref)s | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sed 1d >  %(ref)s.length''' % globals())
    else:
        text = openlenfile('%(ref)s.length' % globals())
        print("\tref.length file already created for %s\n\t\twhich has %s contigs" % (ref, len(text)-1))
    if not op.exists("%(ref)s.length" % globals()):
        print("something went wrong with creating the ref.length file for %s\nexiting %s" % (ref, sys.argv[0]))
        exit()

    # spread contigs across 1000 bed files using the ref.fa.length file
    fcount = make_bedfiles()

    print('\t\tcreated %s bedfiles for %s' % (fcount, ref))


def make_bedfile(lines, fcount, from_orderfile=False):
    """Use ref.fa.length or ref.order file to write contig/chorm, start, stop to bedfile.
    Different than make_bed(): ref.fa.length is 1-based, need to assert zero-based indexing.

    Positional arguments:
    lines - list of tuples - zeroth element is contig/chrom name, first is length of contig/chrom
    fcount - the fcount'th bedfile; used for naming file
    """
    f = get_prereqs(fcount)
    text = []
    with open(f, 'w') as o:
        for line in lines:
            if from_orderfile is False:
                contig, length = line
                text.append("%s\t%s\t%s" % (contig, 0, int(length)-1))  # contig \t start \t stop
            else:
                scaff, start, stop = line
                text.append("%s\t%s\t%s" % (scaff, start, stop))  # scaff \t start \t stop
        o.write("\n".join(text))


def make_bedfiles():
    """Use ref.fa.length file to create bedfiles."""
    text = openlenfile("%s.length" % ref)
    thresh = math.ceil(len(text) / 500)
    lines = []
    fcount = 0
    for count, line in enumerate(text):
        if not line == '':
            lines.append(line.split("\t"))
        if len(lines) == thresh or (count + 1 == len(text)):
            make_bedfile(lines, fcount)
            lines = []
            fcount += 1
    return fcount


def check_beddir():
    """Avoid accidentally using incorrect bedfiles by removing any that exist."""
    bname, beddir = make_beddir()
    files = [f for f in fs(beddir) if f.endswith('.bed')]
    if len(files) > 0:
        text = '\tThere are already existing bedfiles in %s. These will be deleted.' % beddir
        print(Bcolors.WARNING + text + Bcolors.ENDC)
        askforinput()
        for f in files:
            os.remove(f)
        print('\t\tRemoved %s bedfiles.' % len(files))


def main(ref):
    globals().update({'ref': ref})
    
    # warn about overwriting
    check_beddir()

    # get sequence lengths
    find_positions()


if __name__ == "__main__":
    # args
    thisfile, ref = sys.argv
    main(ref)
