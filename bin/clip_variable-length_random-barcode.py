#!/usr/bin/env python2.7

################################################################################
# Demultiplex based on and then clip variable-length barcode region.
# The region is assumed to consist of 3 parts:
# a random N-mer at the beginning, a demultiplexing barcode, and a fixed-sequence spacer.
# The sgRNA sequence follows immediately after this region.
#
# The spacer is used as anchor for finding the demultiplexing barcodes, then entire
# region before the sgRNA gets clipped off.
#
# Input from an unzipped FASTQ pipe.
#
# part of CRISPR / shRNA screen pre-processing pipeline
#
# Kimon Froussios
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2019/05/23
################################################################################

import sys, re, string, csv
from argparse import ArgumentParser, RawDescriptionHelpFormatter


parser = ArgumentParser(description="Variable length trimming of leading n-mer.", formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-o", "--outputdir", type=str, default="./process/fastq", help="Output directory where demultiplexed fastq files will be saved.")
parser.add_argument("-l", "--logFile", type=str, default="./process/fastq/demultiplex_trim.log", help="Report file.")
parser.add_argument("-d", "--demuxOffset", type=int, default=-4, help="Start position of demultiplexing barcode, relative to spacer. Negative integer for upstream of spacer start, positive integer for downstream of spacer end.")
parser.add_argument("-s", "--anchorSeq", type=str, default='TTCCAGCATAGCTCTTAAAC', help="Spacer sequence pattern to anchor.")
parser.add_argument("-b", "--barcodes", type=str, default='barcodes.txt', help="Demultiplexing table, tab-delimited (lane, sample_name, barcode, position). Position is optional, if omitted anchoring will fall back to regex search for the spacer.")
args = parser.parse_args()

# Get barcodes
##############

# Nested dict:  demuxS[lane][barcode] = sample
demuxS = dict() # Sample names
demuxP = dict() # Barcode excpected positions
withPos=False   # Explicit barcode positions available

# Parse barcodes
if args.barcodes:
    with open(args.barcodes, "rt") as bcFile:
        csvreader = csv.DictReader(bcFile, delimiter="\t")
        for i, row in enumerate(csvreader):
            if i == 0:
                if 'position' in row.keys():
                    withPos = True
            if withPos:
                    demuxP[ row['barcode'] ] = row['position']
            m = re.search(r"^([A-Za-z0-9]+_\d+)", row['lane'])
            if not m.group(1) in demuxS.keys():
                demuxS[ m.group(1) ] = dict()
            demuxS[ m.group(1) ][ row['barcode'] ] = row['sample_name']















#
