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

import Levenshtein
from Bio import SeqIO

# Parameters
############

parser = ArgumentParser(description="Demultiplexing with variable length 5' construct of barcode and spacers. Input as SAM format via STDIN. One lane at a time.", formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-L", "--lane", type=str, required=True, help="Name of the lane being piped in, as specified in the barcodes table.")
parser.add_argument("-D", "--outputdir", type=str, default="./process/fastq", help="Output directory where demultiplexed fastq files will be saved.")
parser.add_argument("-l", "--logFile", type=str, default="./process/fastq/demultiplex_trim.log", help="Report file.")
parser.add_argument("-o", "--bcOffset", type=int, default=-4, help="Start position offset of the demultiplexing barcode, relative to the spacer. \
                                                                    Positive for downstream of the spacer end, negative for upstream of the spacer start. (Default -4)")
parser.add_argument("-s", "--anchorSeq", type=str, default='TTCCAGCATAGCTCTTAAAC', help="Spacer sequence to anchor. If using fixed positions specified in the barcodes table, \
                                                                                        this must be a literal sequence, not a regex.")
parser.add_argument("-r", "--anchorRegex", action="store_true", help="--anchorSeq is a regex. (Default False)")
parser.add_argument("-b", "--barcodes", type=str, default='barcodes.txt', help="Demultiplexing table, tab-delimited (lane, sample_name, barcode, position). \
                                                                                Position is 1-based and refers to the start of the anchoring !!SPACER!!, NOT the barcode start! \
                                                                                If omitted, anchoring will fall back to regex search. (default barcodes.txt)")
parser.add_argument("-m", "--bcmm", type=int, default=1, help="Mismatches allowed in matching the demultiplexing barcodes.")
parser.add_argument("-M", "--smm", type=int, default=2, help="Mismatches allowed in matching the spacer sequence.")
args = parser.parse_args()


# Get barcodes
##############

# Dictionaries
demuxS = dict()  # demuxS[barcode] = sample
spacerP = list() # spacerP = [positions]
demuxB = dict()  # demuxB[position] = barcode
withPos = False  # Explicit anchor positions available

# Parse barcodes
if args.barcodes:
    with open(args.barcodes, "rt") as bcFile:
        csvreader = csv.DictReader(bcFile, delimiter="\t")
        for i, row in enumerate(csvreader):
            if i == 0:
                if not ("lane" in row.keys() or "sample_name" in row.keys() or "barcode" in row.keys()):
                    exit("Error: 'lane', 'sample_name', or 'barcode' field is missing from the barcodes table.")
                if 'position' in row.keys():
                    withPos = True              # Spacer start positions have been defined.
            if row['lane'] == args.lane:
                pos = int(row['position']) - 1  # 0-based indexing
                if withPos:
                    if pos < 0:
                        exit(' '.join("Error: Invalid barcode position definition for", row['lane'], row['barcode'], row['sample_name']))
                    if not pos in spacerP:
                        spacerP.append(pos)
                    demuxB[pos] = row['barcode']
                    demuxS[ row['barcode'] ] = row['sample_name']
# Maybe the lane specifications did not match?
if len(spacerP) == 0:
    exit("Error: It looks like no info was parsed from the barcodes table. Does the value of --lane match a value in the table 'lane' column?")

# Clean up lane name
lane = args.lane
if lane[(len(lane)-4):len(lane)] == '.bam':
    lane = lane[0:(len(row['lane'])-4)]         # crop .bam suffix


# Spacer pattern
################
anchor = re.compile(args.anchorSeq) # Pattern matching
anchorLen = len(args.anchorSeq)     # To be used for literal sequence, not pattern

# Open output files
###################
fqOut = dict()
for lane in demuxS.keys():
    for barcode in demuxS[lane].keys():
        laneout = os.path.join(args.outputdir, lane)
        try:
            os.makedirs(laneout)
        except OSError:   # path already exists
            pass
        file = lane + '_' + demuxS[lane][barcode] + '.fastq'
        fqOut[demuxS[lane][barcode]] = open(os.path.join(laneout, file), "w")


# Scan BAM/SAM stream
#####################
for r in sys.stdin:
    #D00689:401:CDM9JANXX:3:1101:1593:1999	4	*	0	0	*	*	0	0	CGGCTNGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAG	BBBBB#<BBFFFFFFFFFFFFFFF<FFFFFFFFFFFFFBFFFFFFFFFFF	B2:Z:NCNNNNCCT	Q2:Z:#<####BBB	BC:Z:GCATTNNNC	RG:Z:CDM9JANXX.3	QT:Z://///###/
    # 0                                     1   2   3   4   5   6   7   8   9                                                   10                                                  11              [12]            [13]            [14]                [15]

    # Split fields
    fields = r.split("\t")
    name = fields[0]
    seq = fields[9]
    qual = fields[10]

    # Anchor
    if withPos:     # With fixed positions
        # Scan through the positions and try to match the anchor.
        for pos in spacerP:
            if hamming(args.anchorSeq, seq[pos:(pos+anchorLen)] <= args.smm:
                
                break

    # without positions

    # Demultiplex

    # Trim

    # Write


    #SeqIO.write(sequences, output_handle, "fastq")

# Close output files
####################
for file in fqOut:
    file.close()



exit(0)
#
