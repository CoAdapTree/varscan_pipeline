"""
### purpose
# filter VariantsToTable output by GQ/globfreq/missing data, keep SNP or INDEL (tipe)
# will also keep biallelic SNPs when REF = N (two ALT alleles)
# for varscan SNP: 
# = filter for non-multiallelic, global MAF >= 1/ploidy_total_across_pops, GQ >= 20, < 25% missing data
# for crisp SNP:
# = filter for non-multiallelic, global MAF >= 1/ploidy_total_across_pops, < 25% missing data
# for INDEL - no filter, just combine (output has multiple rows)
###

### assumes
# gatk VariantsToTable [...] -F TYPE -GF GT -GF GQ [-GF FREQ] --split-multi-allelic
###

### usage
# python filter_VariantsToTable.py SNPorINDEL
# OR
# from filter_VariantsToTable import main as remove_multiallelic
###
"""

import os, sys, pandas as pd, numpy as np, math
from tqdm import tqdm
from coadaptree import uni, pklload
from os import path as op
from collections import Counter


def table(lst):
    c = Counter()
    for x in lst:
        c[x] += 1
    return c


def get_copy(df, cols):
    return df[cols].T.copy()


def get_freq_cutoffs(tablefile):
    pooldir = op.dirname(op.dirname(tablefile))
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)
    poolsamps = pklload(op.join(parentdir, 'poolsamps.pkl'))[pool]
    ploidy = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
    lowfreq = 1/(ploidy * len(poolsamps))
    highfreq = 1 - lowfreq
    return lowfreq, highfreq


def filter_freq(df, tf, tipe, tablefile):
    """
    filter out loci with global MAF < 1/(ploidyPerPop * nPops)
    right now this is unnecessary for varscan when setting pool-level freq to 1/ploidy
    """
    # believe it or not, it's faster to do qual and freq filtering in two steps vs an 'and' statement
    lowfreq, highfreq = get_freq_cutoffs(tablefile)
    print(f'filtering for global frequency ({lowfreq}, {highfreq})...')
    df.reset_index(drop=True, inplace=True)
    freqcols = [col for col in df.columns if '.FREQ' in col]
    copy = get_copy(df, freqcols)
    filtloci = []
    for locus in tqdm(copy.columns):
        freqs = [x for x 
                 in copy[locus].str.rstrip('%').astype('float') 
                 if not math.isnan(x)]  # faster than ...str.rstrip('%').astype('float').dropna()
        if not len(freqs) == 0:
            # avoid loci with all freqs masked (avoid ZeroDivisionError)
            globfreq = sum(freqs)/(100*len(freqs))
            if lowfreq <= globfreq <= highfreq:
                filtloci.append(locus)
    print(f'{tf} has {len(filtloci)} {tipe}s that have global MAF > {lowfreq*100}%')
    df = df[df.index.isin(filtloci)].copy()
    df.index = range(len(df.index))
    return df


def filter_missing_data(df, tf, tipe):
    """
    count np.nan in .FREQ col to assess % missing data
    keep only loci with < 25% missing data
    """
    freqcols = [col for col in df.columns if '.FREQ' in col]
    copy = get_copy(df, freqcols)
    keepers = []
    thresh = math.floor(0.25 * len(freqcols))
    for locus in tqdm(copy.columns):
        # if there is less than 25% missing data:
        # the only time x != x is when x is nan (fastest way to count it)
        count = sum(1 for x in copy[locus] if x != x)
        if count < thresh:
            keepers.append(locus)
    df = df[df.index.isin(keepers)].copy()
    df.index = range(len(df.index))
    return df


def filter_qual(df, tf, tipe, tablefile):
    """mask freqs that have GQ < 20"""
    qualloci = []
    gqcols = [col for col in df.columns if '.GQ' in col]
    thresh = math.ceil(0.75 * len(gqcols))  # assumes len(gqcols) == numpools
    print(f'masking bad freqs for {len(gqcols)} pools...')
    for col in tqdm(gqcols):
        freqcol = col.replace(".GQ", ".FREQ")
        # badloci True if qual < 20
        df.loc[df[col] < 20, freqcol] = np.nan

    print('filtering for missing data ...')
    df = filter_missing_data(df, tf, tipe)

    if len(df.index) > 0:
        print(f'{tf} has {len(df.index)} {tipe}s that have GQ >= 20 and < 25% missing data')
        df = filter_freq(df, tf, tipe, tablefile)
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


def get_refn_snps(df, tipe, ndfs=None):
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
    return (dfs, ndfs)


def recalc_global_freq(df, tf, freqcols):
    """
    For some reason AF reported by crisp is a little off
    - moves crisp AF column to 'crisp_AF'
    - recalulate global AF (alt), save as AF column
    """
    print('Recalculating global freq ... ')
    df.index = range(len(df.index))
    df['crisp_AF'] = df['AF']
    copy = get_copy(df, freqcols)
    for locus in tqdm(copy):
        denom = sum(~copy[locus].isnull())  # num of non-NA
        if denom > 0:
            num = np.nansum(copy[locus])
            res = round(num/denom, 6)
            df.loc[locus,'AF'] = res
        else:
            # all pops have NaN for .FREQ
            # will get filtered later when looking @ missing data
            pass
    return df

def add_freq_cols(df, tf, tipe, tablefile):
    print('Adding in .FREQ columns for crisp file ...')
    # remove bednum from column names so we can pd.concat() later
    bednum = tf.split("file_")[-1].split("_converted")[0]
    df.columns = [col.replace("_" + bednum, "") for col in df.columns]
    # add in a .FREQ column for pool-level freqs
    gtcols = [col for col in df.columns if '.GT' in col]
    print('len(gtcols) = ', len(gtcols))
    freqcols = []
    for col in tqdm(gtcols):
        refcol  = col.replace(".GT", ".REFCOUNT")
        altcol  = col.replace(".GT", ".ALTCOUNT")
        freqcol = col.replace(".GT", ".FREQ")
        freqcols.append(freqcol)
        for alt in uni(df['ALT']):
            df.loc[df['ALT'] == alt, altcol] = df[col].str.count(alt)
        for ref in uni(df['REF']):
            df.loc[df['REF'] == ref, refcol] = df[col].str.count(ref)
        df[freqcol] = df[altcol] / (df[altcol] + df[refcol])
    # remove count cols
    print('Removing unnecessary cols ...')
    df = df[[col for col in df.columns
             if '.REFCOUNT' not in col
             and '.ALTCOUNT' not in col]].copy()
    # recalculate global AF
    df = recalc_global_freq(df, tf, freqcols)
    # sort columns to group data together for each pool
    datacols = sorted([col for col in df.columns if '.' in col])
    othercols = [col for col in df.columns
                 if '.' not in col
                 and col != 'locus'
                 and 'crisp' not in col]
    othercols.insert(othercols.index('AF') + 1, 'crisp_AF')
    df = df[['locus'] + othercols + datacols].copy()
    df.index = range(len(df.index))
    return df


def write_file(tablefile, df, tipe):
    newfile = tablefile.replace(".txt", f"_{tipe}.txt")
    df.to_csv(newfile, index=False, sep='\t')
    print('finished filtering VariantsToTable file: %s' % newfile)


def main(tablefile, tipe, ret=False):
    print('\nstarting filter_VariantsToTable.py for %s' % tablefile)
    tf = op.basename(tablefile)

    # load the data, create a column with CHROM-POS for locusID
    df = pd.read_csv(tablefile, sep='\t')
    print(f'{tf} has {len(df.index)} rows (includes multiallelic)')
    df['locus'] = ["%s-%s" % (contig, pos) for (contig, pos) in zip(df['CHROM'].tolist(), df['POS'].tolist())]

    # determine loci with REF=N but biallelic otherwise
    if tipe == 'SNP':
        dfs, ndfs = get_refn_snps(df, tipe)

        # determine which loci are multiallelic
        loccount = table(df['locus'])
        goodloci = [locus for locus in loccount if loccount[locus] == 1]
        print(f'{tf} has {len(goodloci)} good loci (non-multiallelic)')
    
        # filter df for multiallelic (multiple lines), REF != N
        df = df[df['locus'].isin(goodloci)].copy()
        df = df[df['REF'] != 'N'].copy()

    # filter for tipe, announce num after initial filtering
    df = df[df['TYPE'] == tipe].copy()
    print(f'{tf} has {len(df.index)} good loci of the type {tipe}')
    if len(df.index) == 0:
        if ret is True:
            return df
        else:
            # save
            write_file(tablefile, df, tipe)
    
    # add in loci with REF=N but biallelic otherwise
    if tipe == 'SNP' and len(dfs) > 0:
        print(f'{tf} has {len(ndfs.index)} biallelic {tipe}s with REF=N')
        dfs.append(df)
        df = pd.concat(dfs)
    
    # filter for quality and missing data
    df.index = range(len(df.index))
    if 'varscan' in tf and tipe == 'SNP':
        # if we allow to continue for INDEL, each line is treated as a locus (not true for INDEL)
        df = filter_qual(df, tf, tipe, tablefile)
    if 'crisp' in tf:
        df = add_freq_cols(df, tf, tipe, tablefile)
        print('filtering for missing data ...')
        df = filter_missing_data(df, tf, tipe)
        print(f'{tf} has {len(df.index)} loci with < 25% missing data')
        lowfreq, highfreq = get_freq_cutoffs(tablefile)
        df = df[(df['AF'] <= highfreq) & (df['AF'] >= lowfreq)].copy()
        print(f'{tf} has {len(df.index)} loci with MAF > {lowfreq}')
        df.index = range(len(df.index))
    
    if ret is True:
        return df
    else:
        # save
        write_file(tablefile, df, tipe)


if __name__ == '__main__':
    thisfile, tablefile, tipe = sys.argv

    main(tablefile, tipe)
