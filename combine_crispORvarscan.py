"""
### purpose
# combine output from bedfile output from either CRISP or VarScan
###

### usage
# python combine_crispORvarscan.py pooldir crispORvarscan poolORsamp
###

### assumes
# that all bamfiles were given to samtools in the same order for each bedfile
###
"""

import os, sys, time, random, pandas as pd
from os import path as op
from coadaptree import fs, pklload
from filter_VariantsToTable import main as filtvtt
from start_crispANDvarscan import getfiles, get_bedfiles


def get_varscan_names(df):
    print('renaming varscan columns')
    # get order of samps used to create varscan cmds (same order as datatable)
    pool = op.basename(pooldir)
    samps = pklload(op.join(op.dirname(pooldir), 'poolsamps.pkl'))[pool]
    # create a list of names that varscan gives by default
    generic = ['Sample%s' % (i+1) for i in range(len(samps))]
    # create a map between generic and true samp names
    dic = dict((gen,samp) for (gen,samp) in zip(generic,samps))
    # rename the columns in df
    cols = []
    for col in df:
        if '.' in col:
            gen, rest = col.split(".")
            samp = dic[gen]
            col = '.'.join([samp,rest])
        cols.append(col)
    df.columns = cols
    return df


def checkjobs():
    print('checking jobs')
    parentdir = op.dirname(pooldir)
    pool = op.basename(pooldir)
    ref = pklload(op.join(parentdir, 'poolref.pkl'))[pool]
    samps = fs(op.join(op.dirname(ref),
                       'bedfiles_%s' % op.basename(ref).split(".fa")[0]))
    shdir = op.join(pooldir, 'shfiles/crispANDvarscan')
    # files = {f.sh: f.out, ...}
    files = getfiles(samps, shdir, f"{grep}-{program}")
    return files


def get_types(tablefiles, tipe):
    print(f'starting to filter {len(tablefiles)} tablefiles')
    dfs = [filtvtt(tablefile, tipe, ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)
    
    if program == 'varscan':
        df = get_varscan_names(df)

    filename = op.join(pooldir, f'{program}/{grep}_all_bedfiles_{tipe}.txt')
    df.to_csv(filename, sep='\t', index=False)

    print(f'combined {program} files to {filename}')
    print(f'final {tipe} count = {len(df.index)}')


def get_tables(files):
    print('getting tablefiles')
    tablefiles = [f for f in fs(op.join(pooldir, program))
                  if f.endswith('.txt')
                  and 'all_bedfiles' not in f
                  and 'SNP' not in f
                  and 'INDEL' not in f
                  and grep in f]  # I don't think I need the SNP & INDEL check
    if not len(tablefiles) == len(files):
        print('for some reason tablefiles != files. exiting.')
        exit()
    return tablefiles


def main():
    # make sure all of the crisp jobs have finished
    files = checkjobs()

    # combine table files from output of VariantsToTable
    tablefiles = get_tables(files)
    
    # get SNP and indels
    for tipe in ['SNP', 'INDEL']:
        get_types(tablefiles, tipe)


if __name__ == '__main__':
    # for crisp grep = pool, for varscan grep = pool
    thisfile, pooldir, program, grep = sys.argv 

    main()
