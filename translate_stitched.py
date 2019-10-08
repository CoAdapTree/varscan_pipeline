"""
Translate snps.txt file SNP positions from stitched pos to unstitched pos.

# usage
# python translate_stitched.py snpstable orderfile outfile
#
# if called from command line, args are paths to .txt files
# if imported and main() is called from another script, the snpstable ...
#    and orderfile are each passed as a pandas.DataFrame.ss

# args
# snpsfile = .txt file to be translated
#        - assuming "locus" column (entries for hyphen-separated CHROM-POS)
# ordersfile = ref.order file (no header) with data = (next line)
#        stitched_scaff<tab>unstitched_contig<tab>stitched_start<tab>stitched_stop<tab>contig_length
#        for stitched references, file to translate between ...
#            stitched reference CHROM and unstitched reference CHROM
# outfile = .txt file path to save updated snpstable

# assumes
# .order file has no header
# .order file is of the form ref_scaff<tab>contig_name<tab>start_pos<tab>stop_pos<tab>contig_length
#     positions refer to position of contig within ref_scaff.

# fix
# remove int(pos)-1 from translate() - this was due to the way bedfiles were created.
"""

import sys, pandas as pd
from tqdm import tqdm
from coadaptree import Bcolors


def translate(chrom, pos, order):
    # reduce order to only the contig and range of stitched position
    order = order[(order['stitched_scaff'] == chrom) &
                  (order['stitched_stop'] >= int(pos)-1) &
                  (order['stitched_start'] <= int(pos))
                 ].copy()
    if len(order.index) != 1:
        text = '\tFAIL: did not reduce orderfile correctly'
        print(Bcolors.FAIL + text + Bcolors.ENDC)
        print(chrom, pos)
        print(len(order.index))
        exit()
    order.index = range(len(order.index))
    # get true position
    start,stop = order.loc[0, ['stitched_start', 'stitched_stop']]
    newpos = pos - int(start) + 1
    if not int(start) <= pos <= int(stop):
        text = "\tFAIL: position doesn't make sense"
        print(Bcolors.FAIL + text + Bcolors.ENDC)
        print('\tstart= ', start, 'newpos= ', newpos, 'stop= ', stop)
        exit()
    return order.loc[0, 'unstitched_contig'], newpos


def translate_snps(snps, order):
    """Translate between stiched pos and unstitched position.
    
    # arguments:
    snps = pd.DataFrame
    order = pd.DataFrame
    """
    print('\ttranslating snp positions: stitched -> unstitched')
    new_chroms = []
    new_poss = []
    new_loci = []
    if 'locus' not in snps.columns:
        snps['locus'] = ["%s-%s" % (chrom,pos) for (chrom,pos) in zip(snps['CHROM'],snps['POS'])]
    loci = snps['locus'].tolist()
    for locus in tqdm(loci):
        chrom, pos = locus.split("-")
        new_chrom, new_pos = translate(chrom, int(pos), order.copy())
        new_chroms.append(new_chrom)
        new_poss.append(new_pos)
        new_loci.append("%s-%s" % (new_chrom, new_pos))
    snps['unstitched_chrom'] = new_chroms
    snps['unstitched_pos'] = new_poss
    snps['unstitched_locus'] = new_loci
    
    return snps


def checkfiles(order, snps):
    """Check assumptions for snpstable and orderfile.
    
    # arguments:
    snps = pd.DataFrame read in with chunksize
    order = pd.DataFrame read in without chunksize
    
    # returns orderfile with column names
    """
    print('\tchecking assumptions...')
    if len(order.columns) != 5:
        text = '\tFAIL: there are not 5 columns, as assumed\n\texiting translate_stitched.py'
        print(Bcolors.FAIL + text + Bcolors.ENDC)
        
    order.columns = ['stitched_scaff',
                     'unstitched_contig',
                     'stitched_start',
                     'stitched_stop',
                     'contig_length']
    return order


def main(snps, orderfile, outfile=None):
    text = 'Starting to translate from stitched to unstitched positions'
    print(Bcolors.BOLD + text + Bcolors.ENDC)
    
    # read in order
    order = pd.read_csv(orderfile, sep='\t', header=None, dtype={2:int, 3:int})
    
    # check .order file assumptions
    order = checkfiles(order, snps)
    
    # reduce order to chroms of interest
    order = order[order['stitched_scaff'].isin(snps['CHROM'].tolist())].copy()
    order.index = order['stitched_scaff'].tolist()
    
    # translate snpstable
    translated = translate_snps(snps, order)

    # if called from another script, return the translated dataframe
    if outfile is None:
        return translated
    
    # otherwise, write to table
    print('\twriting to outfile: %s' % outfile)
    translated.to_csv(outfile, sep='\t', index=False)


if __name__ == '__main__':
    thisfile, snpstable, orderfile, outfile = sys.argv
    
    # read in snpstable, assumes snpstable doesn't need to be read in with chunksize
    snps = pd.read_csv(snpstable, sep='\t')

    main(snps, orderfile, outfile)
