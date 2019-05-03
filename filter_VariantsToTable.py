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

import os, sys, pandas as pd, numpy as np, math
from tqdm import tqdm
from coadaptree import uni
from os import path as op
from collections import Counter


def table(lst):
    c = Counter()
    for x in lst:
        c[x] += 1
    return c


def get_copy(df, cols):
    return df[cols].T.copy()


def filter_freq(df, tf, tipe):
    """filter out loci with MAF < 0.05"""
    # believe it or not, it's faster to do qual and freq filtering in two steps vs an 'and' statement
    print('filtering for frequency ...')
    df.reset_index(drop=True, inplace=True)
    freqcols = [col for col in df.columns if '.FREQ' in col]
    copy = get_copy(df, freqcols)
    filtloci = []
    for locus in tqdm(copy.columns):
        freqs = [x for x in copy[locus].str.rstrip('%').astype('float') if not math.isnan(x)]
        globfreq = sum(freqs)/(100*len(freqs))
        if globfreq >= 0.05 and globfreq <= 0.95:
            filtloci.append(locus)
    print(f'{tf} has {len(filtloci)} {tipe}s that have MAF > 5%')
    return df[df.index.isin(filtloci)].copy()


def filter_missing_data(df, tf, tipe):
    """count np.nan in .FREQ col to assess % missing data"""
    freqcols = [col for col in df.columns if '.FREQ' in col]
    copy = get_copy(df, freqcols)
    keepers = []
    thresh = math.floor(0.25 * len(freqcols))
    for locus in tqdm(copy.columns):
        # if there is less than 25% missing data:
        if copy[locus].tolist().count(np.nan) < thresh:
            keepers.append(locus)
    df = df[df.index.isin(keepers)].copy()
    df.index = range(len(df.index))
    return df


def filter_qual(df, tf, tipe):
    """mask freqs that have GQ < 20, keep only loci with < 25% missing data"""
    print('masking bad freqs ...')
    qualloci = []
    gqcols = [col for col in df.columns if '.GQ' in col]
    thresh = math.ceil(0.75 * len(gqcols))  # assumes len(gqcols) == numpools
    for col in tqdm(gqcols):
        freqcol = col.replace(".GQ", ".FREQ")
        badloci = df[col] < 20  # True if qual < 20 else False
        df.loc[badloci, freqcol] = np.nan

    print('filtering for missing data ...')
    df = filter_missing_data(df, tf, tipe)

    if len(df.index) > 0:
        print(f'{tf} has {len(df.index)} {tipe}s that have GQ >= 20 and < 25% missing data')
        df = filter_freq(df, tf, tipe)
        df.index = range(len(df.index))
    else:
        print(f'{tf} did not have any {tipe}s that have GQ >= 20 for >= 75% of pops' +
              '\nnot bothering to filter for freq')
    return df


def adjust_freqs(smalldf):
    """
    for loci with REF=N, set freqs of pools with REF=N in GT to np.nan,
    set alt freqs with respect to the second alt allele
    """
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
    print('starting filter_VariantsToTable.py for %s' % tablefile)
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
    
    # filter for quality and missing data
    df.index = range(len(df.index))
    if 'varscan' in tf and tipe == 'SNP':
        df = filter_qual(df, tf, tipe)

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
