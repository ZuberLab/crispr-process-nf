#!/usr/bin/env sh


## Parameters ##
function usage() {
    echo "Usage:"
    echo "      $0 -i INDIR -l GUIDES_LIB -b DEMUX_TABLE -n COUNTS_OUTDIR -c CONTRASTS_TABLE -m MAGECK_OUTDIR [-1] [-2]"
    echo "      There are many more options and I've lost track. Consult getops in the source code."
		exit 1
}
# Defaults
pad="C"
countthresh=50
spacer='TTCCAGCATAGCTCTTAAAC'
umi=6
demux=4
bcmm=1
bcoffset=-4
smm=2
variable=0
revcomp=0
guideLen=20
FU='/groups/zuber/zubarchive/USERS/Kimon/crispr-process-nf/bin'
PYTHONPATH="$PYTHONPATH:$FU"
# Parse options.
while getopts 'i:l:b:n:r:z:p:s:u:d:O:g:M:A:Vr' flag; do
  case "${flag}" in
    i) indir="${OPTARG}" ;;           # Input directory with unaligned BAMs
    l) library="${OPTARG}" ;;         # sgRNA library
    b) barcodes="${OPTARG}" ;;        # Demultiplexing table: lane \t sample_name \t barcode
    n) countsdir="${OPTARG}" ;;       # Output directory for FASTQs and counts table
    z) countthresh="${OPTARG}" ;;     # Minimum read count per guide (50)
    p) pad="${OPTARG}" ;;             # Padding base ("C")
    s) spacer="${OPTARG}" ;;          # Anchor sequence ('TTCCAGCATAGCTCTTAAAC')
    u) umi="${OPTARG}" ;;             # UMI length (6)
    d) demux="${OPTARG}" ;;           # De-multiplexing barcode length (4)
    O) bcoffset="${OPTARG}" ;;           # De-multiplexing barcode postition relastive to anchor (-4)
    g) guideLen="${OPTARG}" ;;        # Guide length (20). For clipping.
    M) bcmm="${OPTARG}" ;;            # Barcode mismatch allowance (1)
    A) smm="${OPTARG}" ;;            # Anchor mismatch allowance (2)
    V) variable=1 ;;                  # Demultiplex manually, for staggered libraries (no).
    r) revcomp=1 ;;                   # Reverse complement the reads to match the barcodes? (no)
    *) usage ;;
  esac
done

# Check we got the minimum set
if [ -z "$indir" ] || [ -z "$library" ] || [ -z "$barcodes" ] || [ -z "$contrasts" ] || [ -z "$countsdir" ] || [ -z "$mageckdir" ]; then
  echo "-i $indir -l $library -b $barodes -c $contrasts -n $countsdir -m $mageckdir"
  usage
  exit 1
fi

if [ $revcomp -eq 1 ]; then
  revcomp="--reverse_complement"
else
  revcomp=""
fi

## Workflow ##

wait_for_jobs(){
  sleep 60  # seconds, give time to the schedular to put up the task
  sleeptime=120  # ask every 2 mins, for the first 10 mins
  n=1
  while true; do
    if [ $(squeue | grep kimon.fr | grep -c $1) -eq 0 ]; then
      break
    else
      echo sleep another $((sleeptime / 60)) minutes...
      sleep $sleeptime
    fi
    n=$((n + 1))
    if [ $n -eq 5 ]; then
      sleeptime=300  # if still running after 10 mins, ask every 5 mins
    fi
    if [ $n -eq 10 ]; then
      sleeptime=600  # if still running after 30 mins, ask every 10 mins
    fi
  done
}

set -e

libname=$(basename $library)
libname=${libname/.txt/}

if [ $variable -eq 0 ]; then
  echo ''
  echo "Demultiplex BAM, align and count. Using the standard pipeline."
  nextflow run zuberlab/crispr-process-nf $revcomp --inputDir $indir --library $library --padding_base $pad --spacer_seq $spacer --spacer_length ${#spacer} --barcode_demux_length $demux --barcode_random_length $umi --barcode_demux_mismatches $bcmm --barcodes $barcodes --outputDir $countsdir -profile ii2
else
  mkdir -p ${countsdir}/fastq ${countsdir}/fastqc ${countsdir}/counts/${libname} ${countsdir}/aligned/${libname}
  
  echo ''
  echo "Demultiplexing BAM using anchor sequence."
  # Demultiplex
  module load python-levenshtein/0.12.0-foss-2017a-python-2.7.13
  module load pysam/0.14.1-foss-2017a-python-2.7.13
  ${FU}/fileutilities.py T ${indir}/*.bam --loop srun ,--mem=50000 /users/kimon.froussios/crispr-process-nf/bin/demultiplex_by_anchor-pos.py ,-i {abs} ,-D ${countsdir}/fastq ,-l ${countsdir}/fastq/{bas}.log ,-o $bcoffset ,-s $spacer ,-g $guideLen ,-b $barcodes ,-m $bcmm ,-M $smm ,-q 33 ,-Q \&
  module unload python-levenshtein/0.12.0-foss-2017a-python-2.7.13
  module unload pysam/0.14.1-foss-2017a-python-2.7.13
  wait_for_jobs demultip
  
  echo ''
  echo "FastQC (in the background)." # and don't wait for it. I don't need its output for a while.
  module load fastqc/0.11.5-java-1.8.0_121
  ${FU}/fileutilities.py T ${countsdir}/fastq/*/*.fqc --loop srun ,--mem=5000 ,--cpus-per-task 6 fastqc ,-q ,-t 6 ,-f fastq ,-o ${countsdir}/fastqc {abs} \&
  module unload fastqc/0.11.5-java-1.8.0_121
  
  echo ''
  echo "Compressing FASTQ (in the background)."
  ${FU}/fileutilities.py T ${countsdir}/fastq/*/*.fq --loop srun ,--mem=50000 gzip {abs} \&
  
  echo ''
  echo "Guides library to FASTA."
  cw=$(realpath $(pwd))
  cd $(dirname $library)
  srun ${FU}/process_library.R $library C
  cd $cw
  
  echo ''
  echo "Bowtie2 indexing."
  module load bowtie2/2.2.9-foss-2017a
  srun bowtie2-build ${library/.txt/.fasta} ${library/.txt/}
  
  echo ''
  echo "... waiting for gzip to catch up."
  wait_for_jobs gzip
  
  echo ''
  echo "Bowtie2 aligning."
  ${FU}/fileutilities.py T ${countsdir}/fastq/*/*.fq.gz --loop srun ,--mem=10000 ,--cpus-per-task=4 bowtie2 ,-x ${library/.txt/}  ,-U {abs} ,--threads 4 ,-L 20 ,--score-min 'C,0,-1' ,-N 0 ,--seed 42 '2>' ${countsdir}/aligned/${libname}/{bas}.log \> ${countsdir}/aligned/${libname}/{bas}.sam \&
  module unload bowtie2/2.2.9-foss-2017a
  wait_for_jobs bowtie2
  
  echo ''
  echo "Quantifying with featureCounts."
  module load subread/1.5.0-p1-foss-2017a
  ${FU}/fileutilities.py T ${countsdir}/aligned/${libname}/*.sam --loop srun ,--mem=10000 ,--cpus-per-task=4 featureCounts ,-T 4 ,-a ${library/.txt/.saf} ,-F SAF ,-o ${countsdir}/counts/${libname}/{bas}.txt {abs} \&
  module unload subread/1.5.0-p1-foss-2017a
  wait_for_jobs featureC
  
  echo ''
  echo "Combining samples into one table."
  srun --mem=5000 ${FU}/combine_counts.R $library ${countsdir}/counts/${libname}/*.fq.txt > ${countsdir}/counts/${libname}/counts_mageck.txt
  # Fix header. Strip path, strip file extension, strip lane
  mv ${countsdir}/counts/${libname}/counts_mageck.txt ${countsdir}/counts/${libname}/_counts_mageck.txt
  head -n 1 ${countsdir}/counts/${libname}/_counts_mageck.txt | perl -e 'while(<STDIN>){~s/\S+\/(\w{9}_\d_)+//g;~s/\.fq//g;print}' > ${countsdir}/counts/${libname}/counts_mageck.txt
  tail -n +2 ${countsdir}/counts/${libname}/_counts_mageck.txt >> ${countsdir}/counts/${libname}/counts_mageck.txt
  
  echo ''
  echo "MultiQC"
  #wait_for_jobs fastqc  # It should be long finished by now, but better ask.
  module load multiqc/1.3-foss-2017a-python-2.7.13
  srun multiqc -f -x *.run -o ${countsdir}/multiqc ${countsdir}/fastqc ${countsdir}/aligned/${libname} ${countsdir}/counts/${libname}
  module unload multiqc/1.3-foss-2017a-python-2.7.13
  
  echo ''
  echo "Cleaning up intermediate files"
  #rm -r ${countsdir}/fastq
  rm -r ${countsdir}/fastqc
fi

echo ''
echo "Pre-processing finished!"


counts="${countsdir}/counts/${libname}/counts_mageck.txt"
if [[ ! -f $counts ]]; then
  exit 1 "Counts file does not exist."
fi

echo ''
echo "All done!"
