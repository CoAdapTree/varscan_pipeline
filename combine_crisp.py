"""
### purpose
# combine the bedfiles from parallelized crisp runs
###

### usage
# python combine_crisp.py /path/to/pooldir/
###
"""

import os, sys, pandas as pd
from os import path as op
from coadaptree import fs, pklload
from filter_VariantsToTable import main as remove_multiallelic
from start_crisp import getfiles


def checkjobs(pooldir):
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    samps = fs(op.join(op.dirname(ref), 'bedfiles_%s' % op.basename(ref).split(".fa")[0]))
    shdir = op.join(pooldir, 'shfiles/crisp')
    files = getfiles(samps, shdir, 'crisp_bedfile')
    return files


def get_tables(files, pooldir):
    tablefiles = [f for f in fs(op.join(pooldir, 'crisp')) if f.endswith('.txt') and 'all_bedfiles' not in f]
    if not len(tablefiles) == len(files):
        print('for some reason tablefiles != files. exiting.')
        exit()
    dfs = [remove_multiallelic(thisfile, tablefile, ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)
    
    filename = op.join(pooldir, 'crisp/%s_all_bedfiles_snps.txt' % op.basename(pooldir))
    df.to_csv(filename, sep='\t', index=False)
    
    print('combined crisp files to %s' % filename)


def main(pooldir):
    # make sure all of the crisp jobs have finished
    files = checkjobs(pooldir)

    # combine table files from output of VariantsToTable
    get_tables(files, pooldir)


if __name__ == "__main__":
    # args
    thisfile, pooldir = sys.argv

    main(pooldir)
