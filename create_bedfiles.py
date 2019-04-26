"""
### purpose
# create bedfiles from reference so that we can parallelize CRISP
# if an intervals directory exists, use that instead of ref.fa.length file (for stitched refs)
#

### usage
# python create_bedfiles.py /path/to/reference.fasta
###
"""

import sys, os, math
from os import path as op
from coadaptree import fs, makedir


def openlenfile(lenfile):
    with open(lenfile, 'r') as o:
        text = o.read().split("\n")
    return text


def get_prereqs(num):
    bname = op.basename(ref).split(".fa")[0]
    beddir = makedir(op.join(op.dirname(ref), 'bedfiles_%s' % bname))
    f = op.join(beddir, "%s_bedfile_%s.bed" % (bname, str(num).zfill(4)))
    return f


def make_bed(lines, num):
    f = get_prereqs(num)
    with open(f, 'w') as o:
        for contig, start, stop in lines:
            o.write("%s\t%s\t%s\n" % (contig, start, stop))


def make_bed_from_intervals(intdir):
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
    print('created %s bedfiles for %s from interval files' % (len(intfiles), ref))


def make_lenfile():
    refdir = op.dirname(ref)
    intdir = op.join(refdir, 'intervals')
    if op.exists(intdir):
        print('\tusing intervals dir to create bedfiles')
        make_bed_from_intervals(intdir)
        return
    if not op.exists('%(ref)s.length' % globals()):
        print("\tcreating %s.length file (this will take a few minutes)" % op.basename(ref))
        os.system('''cat  %(ref)s | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sed 1d >  %(ref)s.length''' % globals())
    else:
        text = openlenfile('%(ref)s.length' % globals())
        print("\tref.length file already created for %s\n\t\twhich has %s contigs" % (ref, len(text)-1))
    if not op.exists("%(ref)s.length" % globals()):
        print("something went wrong with creating the ref.length file for %s\nexiting %s" % (ref, thisfile))
        exit()

    # spread contigs across 450 bed files
    fcount = make_bedfiles()

    print('created %s bedfiles for %s' % (fcount, ref))


def make_bedfile(lines, fcount):
    f = get_prereqs(fcount)
    with open(f, 'w') as o:
        for line in lines:
            contig, length = line
            o.write("%s\t%s\t%s\n" % (contig, 0, length))  # contig \t start \t stop


def make_bedfiles():
    text = openlenfile("%s.length" % ref)
    thresh = math.ceil(len(text) / 450)
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


def main(thisfile, ref):
    globals().update({'thisfile': thisfile, 'ref': ref})
    # get sequence lengths
    make_lenfile()


if __name__ == "__main__":
    # args
    thisfile, ref = sys.argv
    main(thisfile, ref)
