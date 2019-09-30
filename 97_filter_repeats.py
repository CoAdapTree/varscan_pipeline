"""
Remove SNPs from repeat regions.

# usage
# python 97_filter_repeats.py repeat_file [repeat_dump] [translation_file]
#

# args
# repeat_file = .txt file with header = chrom,start,stop
#               chrom = contig name in unstitched reference
#               start/stop = exact position of repeat region start/stop 
#                     (1-based indexing)
# repeat_dump (optional) = file to save SNPs removed due to position in repeat
# order_file = ref.order file (no header) with data = 
#        stitched_scaff<tab>unstitched_contig<tab>stitched_start<tab>stitched_stop<tab>contig_length
#        for stitched references, file to translate between ...
#            stitched reference CHROM and unstitched reference CHROM
"""

import sys, os
from coadaptree import fs
from os import path as op



def main(repeat, repeat_dump, order):
    
    


if __name__ == '__main__':
    if len(sys.argv) == 2:
        thisfile, repeat_file = sys.argv
        repeat_dump = None
        order = None
    elif len(sys.argv) == 3:
        thisfile, repeat_file, repeat_dump = sys.argv
        order = None
    elif len(sys.argv) == 4:
        thisfile, repeat_file, repeat_dump, order = sys.argv
    
    main(repeat=repeat_file, repeat_dump=repeat_dump, translation=translation)