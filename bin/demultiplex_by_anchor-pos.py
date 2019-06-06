#!/usr/bin/env python2.7

################################################################################
# Demultiplex based on and then clip variable-length barcode region.
#
# The region is assumed to consist of 3 parts:
# a random N-mer at the beginning, a demultiplexing barcode, and a fixed-sequence spacer.
# The sgRNA sequence follows immediately after this region.
#
# The spacer is used as anchor to identify the variable position of the demultiplexing barcodes, 
# then the entire region before the sgRNA gets clipped off.
#
# Input from a BAM file.
# Output to multiple FASTQ files.
#
# part of CRISPR / shRNA screen pre-processing pipeline
#
# Kimon Froussios
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2019/06/05
################################################################################

import sys, os, re, string, csv
from argparse import ArgumentParser, RawDescriptionHelpFormatter

from collections import Counter

import Levenshtein as lev
import pysam


# Parameters
############

parser = ArgumentParser(description="Demultiplexing with variable length 5' construct of barcode and spacers.", formatter_class=RawDescriptionHelpFormatter)
parser.add_argument("-i", "--bam", type=str, required=True, help="Input BAM file. Single-end reads.")
parser.add_argument("-D", "--outputdir", type=str, default="./process/fastq", help="Output directory where demultiplexed fastq files will be saved.")
parser.add_argument("-l", "--tally", type=str, required=False, help="File to write a tally of the reads assigned to each sample. (Default STDOUT)")
parser.add_argument("-o", "--bcOffset", type=int, default=-4, help="Start position offset of the demultiplexing barcode, relative to the spacer. \
                                                                    Positive for downstream of the spacer end, negative for upstream of the spacer start. Negative signs must be excaped. (-4)")
parser.add_argument("-s", "--anchorSeq", type=str, default='TTCCAGCATAGCTCTTAAAC', help="Spacer sequence to anchor. If using fixed positions specified in the barcodes table, \
                                                                                        this must be a literal sequence, not a regex. (TTCCAGCATAGCTCTTAAAC)")
parser.add_argument("-r", "--anchorRegex", action="store_true", help="--anchorSeq is a regex. (False)")
parser.add_argument("-b", "--barcodes", type=str, default='barcodes.txt', help="Demultiplexing table, tab-delimited (lane, sample_name, barcode, position). \
                                                                                Position is 1-based and refers to the start of the anchoring !!SPACER!!, NOT the barcode start! \
                                                                                If omitted, anchoring will fall back to regex search. (./barcodes.txt)")
parser.add_argument("-m", "--bcmm", type=int, default=1, help="Mismatches allowed in matching the demultiplexing barcodes. (1)")
parser.add_argument("-M", "--smm", type=int, default=2, help="Mismatches allowed in matching the spacer sequence. (2)")
parser.add_argument("-q", "--qualOffset", type=int, default=33, help="Base-call quality offset for conversion from pysam to fastq. (33)")
parser.add_argument("-u", "--unmatched", action="store_true", help="Create a FASTQ file for all the reads that did not match the anchor or barcode within the given tolerances. Otherwise they will simply be ignored. (False)")
parser.add_argument("-a", "--abort", type=int, default=30, help="Upper limit for how far into the read to search for the anchor, when no explicit positions are given in the barcodes file. (30)")
args = parser.parse_args()

# Clean up the lane name
lane = os.path.basename(args.bam)
if lane[(len(lane)-4):len(lane)] == '.bam':
    lane = lane[0:(len(lane)-4)]         # crop .bam suffix
    

# Get barcodes
##############

# Dictionaries
demuxS = dict()  # demuxS[barcode] = sample
spacerP = list() # spacerP = [positions]
demuxB = dict()  # demuxB[position] = [barcodes]
withPos = False  # Explicit anchor positions provided

# Parse barcodes
with open(args.barcodes, "rt") as bcFile:
    csvreader = csv.DictReader(bcFile, delimiter="\t")
    for i, row in enumerate(csvreader):
        if i == 0:
            if not ("lane" in row.keys() or "sample_name" in row.keys() or "barcode" in row.keys()):
                raise Exception("Error: 'lane', 'sample_name', or 'barcode' field is missing from the barcodes table.")
            if 'position' in row.keys():
                withPos = True              # Spacer start positions have been defined.
        if row['lane'] == lane or row['lane'] == lane + '.bam':    # Only interested in the details for the lane being demultiplexed by this instance of the script.
            demuxS[ row['barcode'] ] = row['sample_name']
            if withPos:
                pos = int(row['position']) - 1  # 0-based indexing
                if pos < 0:
                    raise ValueError(' '.join("Invalid barcode position definition for", row['lane'], row['barcode'], row['sample_name']))
                if pos not in spacerP:
                    spacerP.append(pos)
                if pos not in demuxB.keys():
                    demuxB[pos] = list()
                demuxB[pos].append(row['barcode'])
            else:
                # Any position is now fair game
                for pos in range(0, args.abort):
                    if pos not in spacerP:
                        spacerP.append(pos)
                    if pos not in demuxB.keys():
                        demuxB[pos] = list()
                    demuxB[pos].append(row['barcode'])

# Maybe the lane specifications did not match?
if len(demuxS) == 0:
    raise Exception("Error: It looks like no info was parsed from the barcodes table. Does the value of --lane match a value in the 'lane' column of the barcodes table?")


# Open output files
###################
fqOut = dict()
laneout = os.path.join(args.outputdir, lane)
for barcode in demuxS.keys():
    try:
        os.makedirs(laneout)
    except OSError:   # path already exists. Hopefully you have permission to write where you want to, so that won't be the cause.
        pass
    file = lane + '_' + demuxS[barcode] + '.fq'
    fqOut[demuxS[barcode]] = open(os.path.join(laneout, file), "w")
unmatched = None
if args.unmatched:
    unmatched = open(os.path.join(args.outputdir, laneout + '_unmatched.fq'), "w")


# Process BAM file
##################
# Spacer pattern
anchor = re.compile(args.anchorSeq) # Pattern matching
anchorLen = len(args.anchorSeq)     # Will be overwritten later if anchorSeq is a regex
# Statistics
counter = Counter()
# Parse SAM
samin = pysam.AlignmentFile(args.bam, "rb", check_sq=False)
for r in samin:
    #D00689:401:CDM9JANXX:3:1101:1593:1999	4	*	0	0	*	*	0	0	CGGCTNGTCAGTATTTTACCAATGACCAAATCAAAGAAATGACTCGCAAG	BBBBB#<BBFFFFFFFFFFFFFFF<FFFFFFFFFFFFFBFFFFFFFFFFF	B2:Z:NCNNNNCCT	Q2:Z:#<####BBB	BC:Z:GCATTNNNC	RG:Z:CDM9JANXX.3	QT:Z://///###/
    # 0                                     1   2   3   4   5   6   7   8   9                                                   10                                                  11              [12]            [13]            [14]                [15]
    counter.update(['total'])
    name = r.query_name
    seq = r.query_sequence
    quals = r.query_qualities
    # Convert qualities to ASCII Phred
    qual = ''
    for q in quals:
        if sys.version_info[0] < 3:
            qual = qual + str(unichr(q + 33))
        else:
            qual = qual + str(chr(q + 33))
    # Find the position of the anchor, within given mismatch tolerance
    anchorFoundAt = None
    for pos in spacerP:     # Scan through predefined positions. T
                            # This also covers the case where no positions were explicitly defined, as all the positions within the allowed range will have been generated instead.
        if args.anchorRegex:    # Just try to match the regex at the required position. Might be less efficient than anchor.search() when all positions are possible. But it's cleaner not creating a separate use case for it.
            m = anchor.match(seq, pos)
            if m:
                anchorFoundAt = m.start()
                anchorLen = m.end() - m.start()
                break
        else:       # Calculate edit distance, not allowing indels in the anchor.
            if lev.hamming(args.anchorSeq, seq[pos:(pos + anchorLen)]) <= args.smm:
                anchorFoundAt = pos
                break
    # Demultiplex, trim
    if anchorFoundAt is not None:   # The anchor could be matched at the given positions with the given mismatch allowance
        anchorEnd = anchorFoundAt + anchorLen
        bcPos = anchorEnd - 1 + args.bcOffset if args.bcOffset > 0 else anchorFoundAt + args.bcOffset
        # Scan through the barcodes expected at this anchor position
        bcfound = False
        for bc in demuxB[anchorFoundAt]:
            bcEnd = bcPos + len(bc)
            if lev.hamming(bc, seq[bcPos:bcEnd]) <= args.bcmm:
                bcFound = True
                trimPos = max(bcEnd, anchorEnd) # Remember, bc can be either up- or down-stream of anchor
                # Print FASTQ entry
                fqOut[demuxS[bc]].write('@' + name + "\n" + seq[trimPos:] + "\n+\n" + qual[trimPos:] + "\n")
                # Keep count
                counter.update(['assigned', demuxS[bc]])
        if not bcFound:
            unmatched.write('@' + name + "\n" + seq + "\n+\n" + qual + "\n")
            counter.update(['BC unmatched'])
    elif args.unmatched:
        unmatched.write('@' + name + "\n" + seq + "\n+\n" + qual + "\n")
        counter.update(['Anchor unmatched'])
samin.close()

# Close output files
####################
for file in fqOut.values():
    file.close()
if args.unmatched:
    unmatched.close()


# Print tally
#############
if args.tally:
    lf = open(args.tally, "w")
    for k,v in counter.most_common():
        lf.write( "\t".join([lane, k, str(v)]) + "\n")
    lf.close()
else:
    for k,v in counter.most_common():
        sys.stdout.write( "\t".join([lane, k, str(v)]) + "\n")


exit(0)
#
