"""
### purpose
# remove multiallelic sites from VariantsToTable output, keep SNP or INDEL (tipe)
# will keep biallelic SNPs when REF = N (two ALT alleles)
# will remove INDELs with REF = N
###

### assumes
# gatk VariantsToTable [...] --split-multi-allelic 
###

### usage
# python filter_VariantsToTable.py SNPorINDEL
# OR
# from filter_VariantsToTable import main as remove_multiallelic
###
"""

import os, sys, pandas as pd, numpy as np
from coadaptree import uni
from os import path as op
from collections import Counter


def table(lst):
    c = Counter()
    for x in lst:
        c[x] += 1
    return c


def adjust_freqs(smalldf):
    gtcols = [col for col in smalldf.columns if 'GT' in col]
    refalt = smalldf.loc[1, 'ALT']
    # adjust freq of first alt with respect to second alt
    for col in gtcols:
        gt = smalldf.loc[1, col]
        if isinstance(gt, str):
            freqcol = col.split(".")[0] + '.FREQ'
            if not gt == 'N/N':
                freq = smalldf.loc[0, freqcol]
                if isinstance(freq, str):
                    if "%" in freq:
                        newfreq = "%s%%" % (100 - float(freq.split("%")[0]))
                        smalldf.loc[0, freqcol] = newfreq
            else:
                # if gt = N/N, adjust to undefined
                smalldf.loc[1, freqcol] = np.nan
        gt2 = smalldf.loc[0, col]
        if isinstance(gt2, str):
            if gt == 'N/N':
                # if gt = N/N, adjust to undefined
                smalldf.loc[0, freqcol] = np.nan
    return smalldf

def main(thisfile, tablefile, tipe, ret=False):
    print('starting %s for %s' % (op.basename(thisfile), tablefile))
    tf = op.basename(tablefile)

    # load the data, create a column with CHROM-POS for locusID
    df = pd.read_csv(tablefile, sep='\t')
    print(f'{tf} has {len(df.index)} rows (includes multiallelic)')
    df['locus'] = ["%s-%s" % (contig, pos) for (contig, pos) in zip(df['CHROM'].tolist(), df['POS'].tolist())]

    # determine loci with REF=N but biallelic otherwise
    if tipe == 'SNP':
        # as far as I can tell, crisp output from convert_pooled_vcf.py will not output REF = N
        ndf = df[df['REF'] == 'N'].copy()
        ndf = ndf[ndf['TYPE'] == tipe].copy()
        ncount = table(ndf['locus'])
        nloci = [locus for locus in ncount if ncount[locus] == 2]
        ndf = ndf[ndf['locus'].isin(nloci)].copy()
        dfs = []
        for locus in uni(ndf['locus']):
            smalldf = ndf[ndf['locus'] == locus].copy()
            if len(smalldf.index) == 2:
                smalldf.index = range(len(smalldf.index))
                smalldf = adjust_freqs(smalldf)
                smalldf.loc[0,'ALT'] = "%s+%s" % (smalldf.loc[0,'ALT'], smalldf.loc[1,"ALT"])
                dfs.append(pd.DataFrame(smalldf.loc[0,:]).T)
        if len(dfs) > 0:
            ndfs = pd.concat(dfs)

    # determine which loci are multiallelic
    loccount = table(df['locus'])
    goodloci = [locus for locus in loccount if loccount[locus] == 1]
    print(f'{tf} has {len(goodloci)} good loci (non-multiallelic)')
    
    # filter df for multiallelic (multiple lines), REF != N, and for SNP
    df = df[df['locus'].isin(goodloci)].copy()
    df = df[df['REF'] != 'N'].copy()
    df = df[df['TYPE'] == tipe].copy()
    print(f'{tf} has {len(df.index)} good loci of the type {tipe}')
    
    # add in loci with REF=N but biallelic otherwise
    if tipe == 'SNP' and len(dfs) > 0:
        print(f'{tf} has {len(ndfs.index)} biallelic {tipe}s with REF=N')
        dfs.append(df)
        df = pd.concat(dfs)

    if ret is True:
        return df
    else:
        # save
        newfile = tablefile.replace(".txt", f"_{tipe}.txt")
        df.to_csv(newfile, index=False, sep='\t')
        print('finished filtering VariantsToTable file: %s' % newfile)


if __name__ == '__main__':
    thisfile, tablefile, tipe = sys.argv

    main(thisfile, tablefile, tipe)
