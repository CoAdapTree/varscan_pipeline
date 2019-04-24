"""
### purpose
# combine output from bedfile output from either CRISP or VarScan
###

### usage
# python combine_crispORvarscan.py pooldir crispORvarscan poolORsamp
###
"""

import os, sys, time, random, pandas as pd
from os import path as op
from coadaptree import fs, pklload
from filter_VariantsToTable import main as remove_multiallelic
from start_crisp import getfiles


def checkjobs():
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    samps = fs(op.join(op.dirname(ref), 'bedfiles_%s' % op.basename(ref).split(".fa")[0]))
    shdir = op.join(pooldir, 'shfiles/%s' % program if program == 'crisp' else 'shfiles/crispANDvarscan')
    files = getfiles(samps, shdir, f"{grep}-{program}")  # files = {f.sh: f.out, ...}
    return files


def get_tables(files):
    tablefiles = [f for f in fs(op.join(pooldir, program))
                  if f.endswith('.txt')
                  and 'all_bedfiles' not in f
                  and grep in f]
    if not len(tablefiles) == len(files):
        print('for some reason tablefiles != files. exiting.')
        exit()
    dfs = [remove_multiallelic(thisfile, tablefile, ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)

    filename = op.join(pooldir, f'{program}/{grep}_all_bedfiles_snps.txt')
    df.to_csv(filename, sep='\t', index=False)

    print(f'combined {program} files to {filename}')
    print(f'final SNP count = {len(df.index)}')


def main():
    # make sure all of the crisp jobs have finished
    files = checkjobs()

    # combine table files from output of VariantsToTable
    get_tables(files)


if __name__ == '__main__':
    # for crisp grep = pool, for varscan grep = pool
    thisfile, pooldir, program, grep = sys.argv 
    
    main()
