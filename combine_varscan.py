"""Filter and combine all gatk VariantsToTable .txt outputs from varscan vcf files.

This file is called only when dependency SLURM_JOBS have completed with exit code 0.
See start_crispANDvarscan.py::create_combine()

### usage
# python combine_crispORvarscan.py pooldir varscan poolORsamp
###

### assumes
# that all bamfiles were given to samtools in the same order for each bedfile
###

### FYI
# filter_VariantsToTable is serially below, as opposed to in parallel after each
# VarScan command, because of the memory requirements. By keeping VarScan .sh files
# at low memory, Priority is affected less and the VarScan jobs are more likely to
# schedule faster than at higher memory requests.
"""

import sys, pandas as pd
from os import path as op
from coadaptree import fs, pklload
from filter_VariantsToTable import main as filtvtt
from start_varscan import getfiles


def checkjobs():
    """
    Make sure previous realigned bamfiles were created without error.
    Avoids unintentionally combining a subset of all final expected files.

    Calls:
    getfiles from start_crispANDvarscan
    """
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


def get_types(tablefiles, tipe, program, pooldir, grep):
    """
    Use filter_VariantsToTable to filter based on tipe {SNP, INDEL}.

    Positional arguments:
    tablefiles - list of paths pointing to the gatk VariantsToTable .txt outputs from varscan or crisp vcf files
    tipe - str; either "SNP" or "INDEL"
    program - str; either "varscan" or "crisp" - used to find and name files
    """
    print(f'starting to filter {len(tablefiles)} tablefiles')
    dfs = [filtvtt(tablefile, tipe, parentdir=op.dirname(pooldir), ret=True)
           for tablefile in tablefiles]
    df = pd.concat(dfs)

    if program == 'varscan':
        df = get_varscan_names(df, pooldir)

    print('writing df to file ...')
    filename = op.join(pooldir, f'{program}/{grep}-{program}_all_bedfiles_{tipe}.txt')
    df.to_csv(filename, sep='\t', index=False)

    print(f'combined {program} files to {filename}')
    print(f'final {tipe} count = {len(df.index)}')


def get_tables(files):
    """Find all existing .txt files, exit if the number doesn't match expectations.

    Positional arguments:
    files - list of shfiles, should be same length as tablefiles (if all jobs are sbatched and done).
    """
    print('getting tablefiles')
    tablefiles = [f for f in fs(op.join(pooldir, program))
                  if f.endswith('.txt')
                  and 'all_bedfiles' not in f
                  and 'SNP' not in f
                  and 'INDEL' not in f
                  and 'REPEATS' not in f
                  and 'PARALOGS' not in f
                  and grep in op.basename(f)]
    if not len(tablefiles) == len(files):
        print('for some reason tablefiles != files. exiting.')
        print('len(tablefiles) = ', len(tablefiles))
        print('len(files) = ', len(files))
        exit()
    return tablefiles


def main():
    # make sure all of the varscan jobs have finished
    files = checkjobs()

    # combine table files from output of VariantsToTable
    tablefiles = get_tables(files)

    # get SNP and indels
    for tipe in ['SNP', 'INDEL']:
        get_types(tablefiles, tipe, program, pooldir, grep)

    # combine repeats and paralogs
    tabledir = op.dirname(tablefiles[0])
    for tipe in ['PARALOGS', 'REPEATS']:
        tablefiles = [f for f in fs(tabledir) if tipe in f and 'all' not in f and f.endswith('.txt')]
        if len(tablefiles) > 0:
            dfs = []
            for t in tablefiles:
                dfs.append(pd.read_csv(t, sep='\t'))
            df = pd.concat(dfs)
            df = get_varscan_names(df, pooldir)
            df.to_csv(op.join(tabledir, f'{op.basename(pooldir)}-{program}_all_bedfiles_{tipe}.txt'),
                      sep='\t', index=False)

if __name__ == '__main__':
    # for crisp grep = pool, for varscan grep = pool
    thisfile, pooldir, program, grep = sys.argv

    main()
