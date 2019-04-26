"""
### purpose
# remove multiallelic sites from VariantsToTable output, keep SNP or INDEL (tipe)
###

### assumes
# gatk VariantsToTable [...] --split-multi-allelic 
###

### usage
# python filter_VariantsToTable.py SNPorINDEL
###
"""

import os, sys
import pandas as pd
from os import path as op
from collections import Counter


def main(thisfile, tablefile, tipe, ret=False):
    print('starting %s for %s' % (op.basename(thisfile), tablefile))
    tf = op.basename(tablefile)

    # load the data, create a column with CHROM-POS for locusID
    df = pd.read_csv(tablefile, sep='\t')
    print(f'{tf} has {len(df.index)} rows (includes multiallelic)')
    df['locus'] = ["%s-%s" % (contig, pos) for (contig, pos) in zip(df['CHROM'].tolist(), df['POS'].tolist())]

    # determine which loci are multiallelic
    loccount = Counter()
    for locus in df['locus']:
        loccount[locus] += 1
    goodloci = [locus for locus in loccount if loccount[locus] == 1]
    print(f'{tf} has {len(goodloci)} good {tipe}s (non-multiallelic)')

    # filter df for multiallelic (multiple lines), and for SNP
    df = df[df['locus'].isin(goodloci)].copy()
    df = df[df['TYPE'] == tipe].copy()
    print(f'{tf} has {len(df.index)} good loci of the type {tipe}')

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
