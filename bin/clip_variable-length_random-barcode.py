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
# Input from an uncompressed SAM stdin stream.
# Output to multiple FASTQ files, one per sample defined in the barcodes table.
#
# part of CRISPR / shRNA screen pre-processing pipeline
#
# Kimon Froussios
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2019/05/23
################################################################################

import sys, os, re, string, csv
from argparse import ArgumentParser, RawDescriptionHelpFormatter

import pysam
from Bio import SeqIO

# Parameters
############

parser = ArgumentParser(description="Variable length trimming of leading n-mer.", formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--bam", type=str, default=sys.stdin, help="Input BAM or SAM filename or stream (default stdin).")
parser.add_argument("-o", "--outputdir", type=str, default="./process/fastq", help="Output directory where demultiplexed fastq files will be saved.")
parser.add_argument("-l", "--logFile", type=str, default="./process/fastq/demultiplex_trim.log", help="Report file.")
parser.add_argument("-d", "--demuxOffset", type=int, default=-4, help="Start position of demultiplexing barcode, relative to spacer. Negative integer for upstream of spacer start, positive integer for downstream of spacer end.")
parser.add_argument("-s", "--anchorSeq", type=str, default='TTCCAGCATAGCTCTTAAAC', help="Spacer sequence pattern to anchor.")
parser.add_argument("-b", "--barcodes", type=str, default='barcodes.txt', help="Demultiplexing table, tab-delimited (lane, sample_name, barcode, position). Position is optional, if omitted anchoring will fall back to regex search for the spacer.")
args = parser.parse_args()

# Get barcodes
##############

# Dictionaries
demuxS = dict() # demuxS[lane][barcode] = sample
demuxP = dict() # demuxP[barcode] = position
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
            lane = row['lane'][0:(len(row['lane'])-4)]         # crop .bam suffix
            if not lane in demuxS.keys():
                demuxS[lane] = dict()
            demuxS[lane][ row['barcode'] ] = row['sample_name']

# Anchor pattern
################
anchor = re.compile(args.anchorSeq)

# Open output files
###################
fqOut = dict()
for lane in demuxS.keys():
    for barcode in demuxS[lane].keys():
        laneout = os.path.join(args.outputdir, lane)
        try:
            os.makedirs(laneout)
        except OSError:
            pass
        file = lane + '_' + demuxS[lane][barcode] + '.fastq'
        fqOut[demuxS[lane][barcode]] = open(os.path.join(laneout, file), "w")



# Connect to BAM/SAM stream
###########################
#bam = pysam.AlignmentFile(args.bam)

# Scan BAM/SAM stream
#####################

#SeqIO.write(sequences, output_handle, "fastq")

# Close output files
####################
#bam.close()



exit(0)
#
