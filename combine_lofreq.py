"""
### purpose
# to combine .txt files from batched LoFreq calls
###

### usage
# python combine_lofreq.py pooldir samp
###
"""

import sys, os, pandas as pd
from os import path as op
from coadaptree import fs
from filter_VariantsToTable import main as remove_multiallelic


def get_tablefiles():
    lodir = op.join(pooldir, 'lofreq')
    return [f for f in fs(lodir) if f.endswith("_table.txt") and samp in f]


def main():
    # shouldn't need to check to see if dependencies have finished, since this job is dependent on afterok

    # get the table files
    tablefiles = get_tablefiles()

    # remove multiallelic and combine
    dfs = [remove_multiallelic(thisfile, tablefile, ret=True) for tablefile in tablefiles]
    df = pd.concat(dfs)

    # write to file
    filename = op.join(pooldir, f'lofreq/{op.basename(pooldir)}-{samp}_all_bedfiles_snps.txt')
    df.to_csv(filename, sep='\t', index=False)
    print('combined crisp files to %s' % filename)


if __name__ == '__main__':
    thisfile, pooldir, samp = sys.argv
    main()