"""
Translate snps.txt file SNP positions from stitched pos to unstitched pos.

# args
# snps = .txt file to be translated
#        - assuming "locus" column (entries for hyphen-separated CHROM-POS)
# order_file = ref.order file (no header) with data = 
#        stitched_scaff<tab>unstitched_contig<tab>stitched_start<tab>stitched_stop<tab>contig_length
#        for stitched references, file to translate between ...
#            stitched reference CHROM and unstitched reference CHROM
# outfile = .txt file path to save updated snpstable
"""

import sys, os, pandas as pd
from os import path as op
from coadaptree import Bcolors, uni


def translate(chrom, pos, order):
    # reduce order to only the contig
    order = order[order['stitched_scaff'] == chrom].copy()
    # isolate stitched position
    order = order[order['stitched_stop'].astype(int) >= int(pos)].copy()
    order = order[order['stitched_start'].astype(int) <= int(pos)].copy()
    if not len(order.index) == 1:
        text = '\tFAIL: did not reduce orderfile correctly'
        print(Bcolors + text + Bcolors.ENDC)
        exit()
    order.index = range(len(order.index))
    # get true position
    start,stop,length = order.loc[0, ['stitched_start', 'stitched_stop', 'contig_length']]
    newpos = int(pos) - int(start) + 1
    if not int(start) >= newpos <= int(stop):
        text = "\tFAIL: position doesn't make sense"
        print(Bcolors + text + Bcolors.ENDC)
        exit()
    return order.loc[0, 'unstitched_contig'], newpos


def translate_snps(snps, order):
    """Translate between stiched pos and unstitched position.
    
    # arguments:
    snps = pd.DataFrame read in with chunksize
    order = pd.DataFrame read in without chunksize
    """
    
    dfs = []
    for chunk in snps:
        new_chroms = []
        new_poss = []
        new_loci = []
        for locus in uni(chunk['locus']):
            chrom, pos = locus.split("-")
            new_chrom, new_pos = translate(chrom, pos)
            new_chroms.append(new_chrom)
            new_poss.append(new_pos)
            new_loci.append("%s-%s" % (new_chrom, new_pos))
        chunk['unstitched_chrom'] = new_chroms
        chunk['unstitched_pos'] = new_poss
        chunk['unstitched_locus'] = new_loci
        dfs.append(chunk)
        # use pandas mode='a' to write_csv without keeping a lot in mem
    
    return pd.concat(dfs)


def checkfiles(snps, order):
    """Check assumptions for snpstable and orderfile.
    
    # arguments:
    snps = pd.DataFrame read in with chunksize
    order = pd.DataFrame read in without chunksize
    
    # returns orderfile with column names
    """
    if len(order.columns) != 5:
        text = '\tFAIL: there are not 5 columns, as assumed\n\texiting translate_stitched.py'
        print(Bcolors.FAIL + text + Bcolors.ENDC)
    if 'locus' not in snps.columns:
        text = '\tFAIL: there should be a column called "locus" in snpstable ...'
        text = text + '\tFAIL: where elements are hyphen-separated CHROM-POS'
        print(Bcolors.FAIL + text + Bcolors.ENDC)
        
    order.columns = ['stitched_scaff',
                     'unstitched_contig',
                     'stitched_start',
                     'stitched_stop',
                     'contig_length']
    return order


def main(snpstable, orderfile, outfile):
    # read in tables
    snps = pd.read_csv(snpstable, sep='\t', chunksize = 10000)
    order = pd.read_csv(order, sep='\t', header=None)

    # check table assumptions
    order = checkfiles()
    
    # translate snpstable
    translated = translate_snps(snps, order)
    
    
    if __name__ != '__main__':
        return outfile


if __name__ == '__main__':
    thisfile, snpstable, orderfile, outfile = sys.argv
    text = 'Starting to translate from stitched to unstitched positions'
    print(Bcolors.BOLD + text + Bcolors.ENDC)
    main(snpstable, orderfile, outfile)