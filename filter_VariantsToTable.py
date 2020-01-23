"""Filter output from gatk VariantsToTable.

### purpose
# filter VariantsToTable output by GQ/globfreq/missing data, keep SNP or INDEL (tipe)
# will also keep biallelic SNPs when REF = N (two ALT alleles) that pass filters
# for varscan SNP:
# = filter for non-multiallelic, global MAF >= 1/ploidy_total_across_pops, GQ >= 20, < 25% missing data

###
# if number of samples == 1, then no SNPs will be filtered for MAF

### assumes
# gatk VariantsToTable [...] -F TYPE -GF GT -GF GQ -GF FREQ --split-multi-allelic

### usage
## to filter VariantsToTable output file:
# python filter_VariantsToTable.py VariantsToTable_output.txt SNP
# python filter_VariantsToTable.py VariantsToTable_output.txt INDEL
## to remove paralogs and/or remove repeats and/or translate stitched positions:
## (parentdir is used to find .pkl files to determine if apporpriate for pool)
# python filter_VariantsToTable.py VariantsToTable_output.txt SNP parentdir
## OR within another module
# from filter_VariantsToTable import main as remove_multiallelic

### fix
# if a tablefile has zero rows, it will not output a PARALOGS or REPEATS file
"""

import sys, pandas as pd, numpy as np, math
import translate_stitched
from tqdm import tqdm
from coadaptree import uni, pklload
from os import path as op
from collections import Counter


def table(lst):
    """
    Count each item in a list.
    
    Returns:
    - Counter() with key = item, val = count
    """
    c = Counter()
    for x in lst:
        c[x] += 1
    return c


def get_copy(df, cols):
    """
    Transpose dataframe using specific columns (that will be index after transformation).
    Doing so helps speed things up.
    """
    return df[cols].T.copy()


def get_freq_cutoffs(tablefile):
    """
    Determine MAF using ploidy and the number of samples per pool.
    Sums across ploidy values for a given pool to determin MAF.
    
    Assumes:
    - equal ploidy across samples/pools
    
    Positional arguments:
    tablefile - path to VariantsToTable output - used to find ploidy etc
    
    Returns:
    lowfreq - minimum allele freq to keep (MAF)
    highfreq - maximum allele freq to keep (1-MAF)
    """
    pooldir = op.dirname(op.dirname(tablefile))
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)
    poolsamps = pklload(op.join(parentdir, 'poolsamps.pkl'))[pool]
    ploidy = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
    lowfreq = 1/sum(ploidy.values())
    if lowfreq == 1 or len(poolsamps) == 1:
        # for megagametophyte data
        lowfreq = 0
    pklfile = op.join(parentdir, 'maf.pkl')
    if op.exists(pklfile):
        lowfreq = float(pklload(pklfile))
    highfreq = 1 - lowfreq
    return lowfreq, highfreq


def filter_freq(df, tf, tipe, tablefile):
    """
    Filter out loci with global MAF < 1/(total_ploidy_across_pools).
    Right now this is unnecessary for varscan when setting pool-level freq to 1/ploidy.
    
    Positional arguments:
    df - pandas.dataframe; VariantsToTable output
    tablefile - path to VariantsToTable output - used to find ploidy etc
    tf - str; basename of tablefile
    tipe - str; one of either "SNP" or "INDEL"
    
    Returns:
    df - pandas.dataframe; freq-filtered VariantsToTable output
    """
    # believe it or not, it's faster to do qual and freq filtering in two steps vs an 'and' statement
    lowfreq, highfreq = get_freq_cutoffs(tablefile)
    print(f'filtering for global frequency ({lowfreq}, {highfreq})...')
    df.reset_index(drop=True, inplace=True)
    
    # prep for filtering
    freqcols = [col for col in df.columns if '.FREQ' in col]
    pool = op.basename(op.dirname(op.dirname(tablefile)))
    parentdir = op.dirname(op.dirname(op.dirname(tablefile)))
    ploidy = pklload(op.join(parentdir, 'ploidy.pkl'))[pool]
    
    # carry on with poolseq datas
    filtloci = []
    afs = []
    copy = get_copy(df, freqcols)
    for locus in tqdm(copy.columns):
        freqs = dict((samp.replace(".FREQ",""),freq) for (samp,freq)
                     in copy[locus].str.rstrip('%').astype('float').items()
                     if not math.isnan(freq))  # faster than .str.rstrip('%').astype('float').dropna()
        if len(freqs) > 0:  # avoid loci with all freqs masked (avoid ZeroDivisionError)
            # calc globfreq using the samps/ploidy that are present for this locus
            globfreq = sum([ploidy[samp]*(freq/100)
                            for (samp,freq) in freqs.items()]) / sum([ploidy[samp] for samp in freqs])
            if lowfreq <= globfreq <= highfreq:
                filtloci.append(locus)
                # since we're going in order of rows in df ...
                # ... we can use afs to replace AF col later since we reduce df to filtloci
                afs.append(globfreq)
                # which is about 40x faster than: df.loc[locus, 'AF'] = globfreq
    print(f'{tf} has {len(filtloci)} {tipe}s that have global MAF > {lowfreq*100}%')
    df = df[df.index.isin(filtloci)].copy()
    df.index = range(len(df.index))
    df['AF'] = afs
    return df


def filter_missing_data(df, tf, tipe):
    """
    Remove loci with < 25% missing data.
    Count np.nan in .FREQ col to assess % missing data.
    
    Positional arguments:
    df - pandas.dataframe; VariantsToTable output
    tf - str; basename of tablefile
    tipe - str; one of either "SNP" or "INDEL"
    
    Returns:
    df - pandas.dataframe; missing data-filtered VariantsToTable output
    """
    freqcols = [col for col in df.columns if '.FREQ' in col]
    copy = get_copy(df, freqcols)
    keepers = []
    # else statement for running single pop (megagamtophyte) through:
    thresh = math.floor(0.25 * len(freqcols)) if len(freqcols) > 1 else 1
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
    """
    mask freqs that have GQ < 20.
    
    Positional arguments:
    df - pandas.dataframe; VariantsToTable output
    tf - str; basename of tablefile
    tipe - str; one of either "SNP" or "INDEL"
    
    Returns: pandas.dataframe; quality-filtered VariantsToTable output
    - FREQ and GT are masked (np.nan) if GQ < 20
    """
    gqcols = [col for col in df.columns if '.GQ' in col]
    print(f'masking bad freqs for {len(gqcols)} pools...')
    for col in tqdm(gqcols):
        freqcol = col.replace(".GQ", ".FREQ")
#         gtcol = col.replace(".GQ", ".GT")  # pretty sure this is depricated
        # badloci True if qual < 20
#         df.loc[df[col] < 20, [freqcol, gtcol]] = np.nan
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
#         df = drop_freq_cols(df)
    return df


def adjust_freqs(smalldf):
    """
    For loci with REF=N, set freqs of pools with REF=N in GT to np.nan.
    Set alt freqs with respect to the second alt allele.
    
    Positional arguments:
    smalldf - pandas.dataframe; df with only REF=N
    
    Returns:
    ndf - smalldf with adjusted freqs in zeroth row
    """
    gtcols = [col for col in smalldf.columns if 'GT' in col]

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
    """
    Isolate polymorphisms with REF=N but two ALT single nuleodite alleles.
    
    Positional arguments:
    df - pandas.dataframe; current filtered VariantsToTable output
    
    Returns:
    dfs - list of loci (pandas.dataframes) with REF=N and two ALT alleles, counts with respect to second ALT
    ndfs - return from pd.conat(dfs)
    """
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


def write_file(tablefile, df, tipe):
    """Write filtered pandas.dataframe to file using args to create file name."""
    newfile = tablefile.replace(".txt", f"_{tipe}.txt")
    df.to_csv(newfile, index=False, sep='\t')
    print('finished filtering VariantsToTable file: %s' % newfile)


def get_varscan_names(df, pooldir):
    """Convert generic sample/pool names from varscan to something meaningful."""
    print('renaming varscan columns ...')
    # get order of samps used to create varscan cmds (same order as datatable)
    pool = op.basename(pooldir)
    print('pklfile = ', op.join(op.dirname(pooldir), 'poolsamps.pkl'))
    samps = pklload(op.join(op.dirname(pooldir), 'poolsamps.pkl'))[pool]
    print('len(samps) for %s = ' % pool, len(samps))
    # create a list of names that varscan gives by default
    generic = ['Sample%s' % (i+1) for i in range(len(samps))]
    # create a map between generic and true samp names
    dic = dict((gen, samp) for (gen, samp) in zip(generic, samps))
    # rename the columns in df
    cols = []
    for col in df:
        if '.' in col:
            gen, rest = col.split(".")
            samp = dic[gen]
            col = '.'.join([samp, rest])
        cols.append(col)
    df.columns = cols
    return df


def load_data(tablefile):
    """
    Load the VariantsToTable output.
    
    Positional arguments:
    tablefile - path to VariantsToTable output - used to find ploidy etc
    
    Returns:
    df - pandas.dataframe; VariantsToTable output
    tf - basename of tablefile
    """
    tf = op.basename(tablefile)
    pooldir = op.dirname(op.dirname(tablefile))

    # load the data, create a column with CHROM-POS for locusID
    df = pd.read_csv(tablefile, sep='\t')
    print(f'{tf} has {len(df.index)} rows (includes multiallelic)')
    df['locus'] = ["%s-%s" % (contig, pos) for (contig, pos) in zip(df['CHROM'].tolist(), df['POS'].tolist())]
    df = get_varscan_names(df, pooldir)

    return df, tf, pooldir


def keep_snps(df, tf):
    """
    Count CHROM-POS (locus) and keep only those with one ALT.
    
    Positional arguments:
    df - pandas.dataframe; currently filtered VariantsToTable output
    tf - basename of path to VariantsToTable output

    Returns:
    df - pandas.dataframe; non-multiallelic-filtered VariantsToTable output
    """
    loccount = table(df['locus'])
    goodloci = [locus for locus in loccount if loccount[locus] == 1]
    print(f'{tf} has {len(goodloci)} good loci (non-multiallelic)')

    # filter df for multiallelic (multiple lines), REF != N
    df = df[df['locus'].isin(goodloci)].copy()
    df = df[df['REF'] != 'N'].copy()
    return df


def filter_type(df, tf, tipe):
    """Keep onli loci called a SNP by program."""
    df = df[df['TYPE'] == tipe].copy()
    print(f'{tf} has {len(df.index)} good loci of the type {tipe}')
    return df


def remove_paralogs(snps, parentdir, snpspath, pool):
    """
    Remove sites from snptable that are thought to have multiple gene copies align to this position.
    
    # assumes
    # paralog file has 'CHROM' and 'locus' in the header (best if this is the only data, reads in quicker)
    #   where CHROM is the reference chromosome/scaffold
    #   where locus is hyphen-separated CHROM-POS
    
    # paralog file is created from calling SNPs on haplotype data as diploid
    #   no need to worry about translating stiched -> unstitched if SNPs called on same reference.
    """
    parpkl = op.join(parentdir, 'paralog_snps.pkl')
    if op.exists(parpkl):
        # read in paralogfile
        paralogdict = pklload(parpkl)
        if paralogdict[pool] is not None:
            print('Removing paralogs sites ...')
            paralogs = pd.read_csv(paralogdict[pool], sep='\t')
            # remove and isolate paralogs from snps
            truths = snps['locus'].isin(paralogs['locus'])
            found_paralogs = snps[truths].copy()
            snps = snps[~truths].copy()
            snps.index = range(len(snps.index))

            # write paralogs to a file
            parafile = snpspath.replace(".txt", "_PARALOGS.txt")
            found_paralogs.to_csv(parafile, sep='\t', index=False)
            print(f'{op.basename(snpspath)} has {len(snps.index)} non-paralog SNPs')
    return snps


def remove_repeats(snps, parentdir, snpspath, pool):
    """
    Remove SNPs that are found to be in repeat-masked regions.
    
    # assumes
    # that the positions have been translated BEFORE removing repeats
        # took forever to create unstitched repeat regions, don't want to translate repeat file
        # this way I can just use unstitched chrom if reference is stitched
    # repeat file has a header ('CHROM', 'start', 'stop')
    # start and stop positions of repeat regions are 1-based
    """
    reppkl = op.join(parentdir, 'repeat_regions.pkl')
    if op.exists(reppkl):
        # read in repeat regions
        repeatdict = pklload(reppkl)
        if repeatdict[pool] is not None:
            print('Removing repeat regions ...')
            # if user selected translation be applied to this pool
            repeats = pd.read_csv(repeatdict[pool], sep='\t')
            # figure out if data is from stitched or not
            if 'unstitched_chrom' in snps.columns:
                # then the snps have been translated: stitched -> unstitched
                chromcol = 'unstitched_chrom'
                poscol = 'unstitched_pos'
                print('\tsnps have been translated')
            else:
                # otherwise SNPs were called on unstitched reference
                chromcol = 'CHROM'
                poscol = 'POS'
                print('\tsnps have not been translated')
            # reduce repeats to the chroms that matter (helps speed up lookups)
            repeats = repeats[repeats['CHROM'].isin(snps[chromcol].tolist())].copy()

            # isolate SNPs in repeat regions
            repeat_snps = []
            for chrom in tqdm(uni(snps[chromcol])):
                reps = repeats[repeats['CHROM'] == chrom].copy()
                mysnps = snps[snps[chromcol] == chrom].copy()
                if len(reps.index) > 0 and len(mysnps.index) > 0:
                    for row in mysnps.index:
                        pos = snps.loc[row, poscol]  # index is maintained from snps to mysnsps
                        df = reps[reps['stop'].astype(int) >= int(pos)].copy()
                        df = df[df['start'].astype(int) <= int(pos)].copy()
                        if len(df.index) > 0:
                            assert len(df.index) == 1
                            repeat_snps.append(row)

            # save repeats
            print(f'\tSaving {len(repeat_snps)} repeat regions')
            repeat_path = snpspath.replace(".txt", "_REPEATS.txt")
            myrepeats = snps[snps.index.isin(repeat_snps)].copy()
            myrepeats.to_csv(repeat_path, sep='\t', index=False)

            # remove SNPs in repeat regions
            snps = snps[~snps.index.isin(repeat_snps)].copy()
            snps.index = range(len(snps.index))

            print(f'{op.basename(snpspath)} has {len(snps.index)} SNPs outside of repeat regions')

    return snps

def translate_stitched_to_unstitched(df, parentdir, pool):
    """See if user asked regions to be translated from stitched genome to unstitched.

    # assumes
    # that this is run BEFORE removing repeats
    """
    orderpkl = op.join(parentdir, 'orderfile.pkl')
    if op.exists(orderpkl):
        orderdict = pklload(orderpkl)
        if oderdict[pool] is not None:
            # if user selected translation be applied to this pool
            orderfile = orderdict[pool]
            df = translate_stitched.main(df.copy(), orderfile)
    return df


def main(tablefile, tipe, parentdir=None, ret=False):
    print('\nstarting filter_VariantsToTable.py for %s' % tablefile)

    # load the data
    df, tf, pooldir = load_data(tablefile)

    # determine loci with REF=N but biallelic otherwise
    if tipe == 'SNP':
        dfs, ndfs = get_refn_snps(df, tipe)

        # determine which loci are multiallelic
        df = keep_snps(df, tf)

    # filter for tipe, announce num after initial filtering
    df = filter_type(df, tf, tipe)

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

    # look for filtering options called at 00_start.py
    if parentdir is not None and tipe == 'SNP':
        # translate stitched (if called at 00_start)
        # translate before filtering so that REPEAT and PARALOG files are translated
        # (takes longer than if translating after filtering, obv)
        df = translate_stitched_to_unstitched(df.copy(),
                                              parentdir,
                                              op.basename(pooldir))  # not coded to translate INDELs

        # remove repeats (if called at 00_start) - want to remove repeats before paralogs
        df = remove_repeats(df.copy(),
                            parentdir,
                            tablefile,
                            op.basename(pooldir))

        # remove paralog SNPs (if called at 00_start)
        df = remove_paralogs(df.copy(), parentdir, tablefile, op.basename(pooldir))

    if ret is True:
        return df
    else:
        # save
        write_file(tablefile, df, tipe)


if __name__ == '__main__':
    if len(sys.argv) == 3:
        thisfile, tablefile, tipe = sys.argv
        parentdir=None
    elif len(sys.argv) == 4:
        # use parentdir to store pkl files with paths to translate stitched, repeat regions, and paralogs
        thisfile, tablefile, tipe, parentdir = sys.argv

    main(tablefile, tipe, parentdir)
