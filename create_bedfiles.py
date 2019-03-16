### purpose
# create bedfiles from reference so that we can parallelize CRISP
#

### usage
# python create_bedfiles.py /path/to/reference.fasta
###

### imports
import sys, os, math
from os import path as op
###

def makedir(DIR):
    if not op.exists(DIR):
        os.makedirs(DIR)
    return DIR

def openlenfile(lenfile):
    with open(lenfile,'r') as o:
        text = o.read().split("\n")
    return text

def make_lenfile():
    if not op.exists('%(ref)s.length' % globals()):
        os.system('echo "creating %s.length file (this will take a few minutes)"' % op.basename(ref))
        os.system('''cat  %(ref)s | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' | sed 1d >  %(ref)s.length''' % globals())
    else:
        text = openlenfile('%(ref)s.length' % globals())
        os.system('echo -e "ref.length file already created for %s\n\twhich has %s contigs"' % (ref,len(text)-1))
    if not op.exists("%(ref)s.length" % globals()):
        os.system('echo -e "something went wrong with creating the ref.length file for %s\nexiting %s"' % (ref,thisfile))
        exit()
    
def make_bedfile(lines,fcount):
    bname = op.basename(ref).replace(".fasta","")
    beddir = makedir( op.join(op.dirname(ref),'bedfiles_%s' % bname) )
    f = op.join(beddir,"%s_bedfile_%s.bed" % (bname,str(fcount).zfill(2)))
    with open(f,'w') as o:
        for line in lines:
            contig,length = line
            o.write("%s\t%s\t%s\n" % (contig,1,length)) # contig \t start \t stop
    
def make_bedfiles():
    text = openlenfile("%s.length" % ref)
    thresh = math.ceil(len(text) / 15)
    lines = []
    fcount = 0
    for count,line in enumerate(text):
        if not line == '':
            lines.append(line.split("\t"))
        if len(lines) == thresh or (count + 1 == len(text)):
            make_bedfile(lines,fcount)
            lines = []
            fcount += 1
    return fcount

def main(thisfile,ref):
    globals().update({'thisfile':thisfile,'ref':ref})
    # get sequence lengths
    make_lenfile()
    
    # spread contigs across 15 bed files
    fcount = make_bedfiles()
    
    os.system('echo created %s bedfiles for %s' % (fcount,ref))

if __name__ == "__main__":
    # args
    thisfile, ref = sys.argv
    main(thisfile,ref)