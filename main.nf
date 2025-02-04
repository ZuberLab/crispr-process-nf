#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ================================================================
     crispr-process-nf
    ================================================================
     DESCRIPTION

     Process CRISPR and shRNA functional genetic screening data.

     Usage:
     nextflow run zuberlab/crispr-process-nf

     Options:
        --inputDir                  Input directory containing raw files.
                                    (default: '01_raw')

        --outputDir                 Output directory for processed files.
                                    (default: 'results')

        --library                   Path to sgRNA / shRNA library file.
                                    (default: 'library.txt')
                                    The following columns are required:
                                        - id:       unique name of sgRNA / shRNA
                                        - group:     gene targeted by sgRNA / shRNA
                                        - sequence: nucleotide sequence of sgRNA / shRNA

         --barcodes                 Path to file containing barcodes for demultiplexing.
                                    (default: 'barcodes.txt')
                                    The following columns are required:
                                        - lane:         name of BAM / FASTQ input file
                                        - sample_name:  name of demultiplexed sample 
                                        - barcode:      sequence of the sample barcode, 
                                                        must be unique for each sample


        --barcode_demux_mismatches  Number of mismatches allowed during demultiplexing
                                    of barcode. (default: 1)

        --random_barcode_length     Number of nucleotides in random barcode.
                                    (default: 2)

        --spacer_seq                Nucleotide sequence in spacer sequence between
                                    barcodes and sgRNA / shRNA sequence. 
                                    (default: TTCCAGCATAGCTCTTAAAC)

        --max_guide_length          Number of nucleotides in guide sequence. (default: 21)

        --rc_guide_seq               Reverse complement guide sequence. (default: false)

        --padding                   Nucleotides used for padding if sgRNA / shRNA are of
                                    unequal length. Corresponds to nucleotides downstream of 
                                    the sgRNA (post guide sequence). Must be one of G, C, T, 
                                    and A. (default: GTT)

        --add_unknown_to_fastqc     Add unknown sequences (no barcode match during demultiplexing) 
                                    to multiQC report. (default: false)

     Profiles:
        standard                    local execution
        singularity                 local execution with singularity
        cbe                         SLURM execution with singularity on CBE cluster

     Docker:
     zuberlab/crispr-nf:0.6

     Author:
     Florian Andersch (florian.andersch@imp.ac.at)
     based on previous work by: Jesse J. Lipp & Tobias Neumann
    """.stripIndent()
}

params.help = false

if (params.help) {
    helpMessage()
    exit 0
}

log.info ""
log.info " parameters "
log.info " ======================"
log.info " input directory                      : ${params.inputDir}"
log.info " output directory                     : ${params.outputDir}"
log.info " library file                         : ${params.library}"
log.info " barcode file                         : ${params.barcodes}"
log.info " random barcode (nt)                  : ${params.random_barcode_length}"
log.info " spacer (nt)                          : ${params.spacer_seq}"
log.info " demultiplex mismatches               : ${params.barcode_demux_mismatches}"
log.info " guide padding bases                  : ${params.padding}"
log.info " max guide length                     : ${params.max_guide_length}"
log.info " reverse complement guide sequence    : ${params.rc_guide_seq}"
log.info " add unknown to FastQC                : ${params.add_unknown_to_fastqc}"
log.info " ======================"
log.info ""

Channel
    .fromPath( "${params.inputDir}/*.bam" )
    .map { file -> tuple( file.baseName, file ) }
    .set { rawBamFiles }

Channel
    .fromPath( "${params.inputDir}/*.{fastq,fq}.gz" )
    .map { file -> tuple( file.baseName.replaceAll(/\.fastq|\.fq/, ''), file ) }
    .set { rawFastqFiles }

Channel
    .fromPath(params.barcodes)
    .set { processBarcodeFiles }

Channel
    .fromPath(params.barcodes)
    .set { barcodeFile }

Channel
    .fromPath(params.library)
    .into { processLibraryFile ; combineLibraryFile }

process bam_to_fastq {

    tag { lane }

    input:
    set val(lane), file(bam) from rawBamFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into fastqFilesFromBam

    script:
    """
    samtools fastq -@ ${task.cpus} ${bam} -2 read_pairs_not_used.fastq > ${lane}.fastq
    rm read_pairs_not_used.fastq
    pigz -p ${task.cpus} ${lane}.fastq
    """
}

fastqFilesFromBam
    .mix(rawFastqFiles)
    .set { fastqFiles }


process trim_random_barcode {

    tag { lane }

    input:
    set val(lane), file(fastq) from fastqFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into randomBarcodeTrimmedFiles

    script:
    """
    mv ${fastq} input.fastq.gz
    spacer="${params.spacer_seq}"
    length_spacer=\${#spacer}
    cutadapt -u ${params.random_barcode_length} --overlap \${length_spacer} -g ${params.spacer_seq} --action=none -j ${task.cpus} input.fastq.gz -o ${lane}.fastq.gz
    """
}

process process_barcodes {

    tag { barcodes.baseName }

    input:
    file(barcodes) from processBarcodeFiles

    output:
    file("*.fasta") into demuxBarcodeFiles

    script:
    """
    process_barcodes.R ${barcodes}
    """
}

demuxBarcodeFiles
    .flatten()
    .map { file -> [file.baseName, file] }
    .concat(randomBarcodeTrimmedFiles)
    .groupTuple()
    .set { demuxFiles }

process demultiplex {

    tag { lane }

    input:
    set val(lane), file(files) from demuxFiles

    output:
    set val(lane), file('*.fq.gz') into splitFiles

    script:
    """
    cutadapt -j ${task.cpus} -e ${params.barcode_demux_mismatches} --no-indels -g file:${files[0]} --action=none -o "${lane}#{name}.fq.gz" ${files[1]}
    """
}

def ungroupTuple = {
    def result = []
    def name = it[0]
    it[1].each { result << [name, it] }
    return result
 }

splitFiles
    .flatMap { it -> ungroupTuple(it) }
    .map { lane, file -> tuple(lane, file.name.replaceAll(/\.fq\.gz/, ''), file) }
    .into { flattenedSplitFilesUnknown; fastqcSplitFilesUnknown }

flattenedSplitFilesUnknown
    .filter { it =~ /^(?!.*unknown).*$/ }
    .into { flattenedSplitFiles; fastqcSplitFiles }

process trim_barcode_and_spacer {

    tag { id }

    publishDir path: "${params.outputDir}/fastq/${lane}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(lane), val(id), file(fastq) from flattenedSplitFiles
    each file(barcodes) from barcodeFile

    output:
    set val(lane), val(id), file("${id}.fastq.gz") into spacerTrimmedFiles

    script:
    sample = id.replaceAll(/^.*?#/, '')
    """
    barcode=\$(awk -F'\\t' '\$2 == "${sample}" {print \$3}' ${barcodes})
    barcode_spacer="\${barcode}${params.spacer_seq}"
    length_barcode_spacer=\${#barcode_spacer}

    mv ${fastq} input.fastq.gz
    cutadapt input.fastq.gz -j ${task.cpus} -u \${length_barcode_spacer} -o ${id}.fastq.gz -l ${params.max_guide_length}
    fastqc -q ${id}.fastq.gz
    """
}

process process_library {

    tag { library.baseName }

    input:
    file(library) from processLibraryFile

    output:
    file("${library.baseName}.saf") into librarySafFile
    file("${library.baseName}.fasta") into libraryFastaFile

    script:
    """
    process_library.R ${library} ${params.padding}
    """
}

process bowtie_index {

    tag { library_fasta.baseName }

    input:
    file(library_fasta) from libraryFastaFile

    output:
    file("bt2") into bt2Index

    script:
    """
    mkdir -p bt2
    bowtie2-build ${library_fasta} bt2/index
    """
}

process align {

    tag { id }

    input:
    set val(lane), val(id), file(fastq) from spacerTrimmedFiles
    each file(index) from bt2Index

    output:
    set val(lane), val(id), file("${id}.sam") into alignedFiles
    file "${id}.log" into alignResults

    script:
    if (params.rc_guide_seq) {
        rc = "--nofw"
    } else {
        rc = "--norc"
    }

    """
    bowtie2 \
        --threads \$((${task.cpus})) \
        -x ${index}/index \
        -L ${params.max_guide_length} \
        --score-min 'C,0,-1' \
        -N 0 \
        --seed 42 \
        $rc \
        <(zcat ${fastq}) 2> ${id}.log > ${id}.sam
        
    """
}

alignedFiles
    .map { lane, id, file -> tuple(lane, file) }
    .groupTuple()
    .set { groupedAlignedFiles }

process count {

    tag { lane }

    publishDir path: "${params.outputDir}/counts/${saf.baseName}/${lane}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(lane), file(sams) from groupedAlignedFiles
    each file(saf) from librarySafFile

    output:
    file("${lane}.txt") into countedFiles
    file("${lane}.txt.summary") into featureCountsResults

    script:
    """
    featureCounts \
        -T ${task.cpus} \
        -a ${saf} \
        -F SAF \
        -o ${lane}.txt \
        ${sams}
    """
}

process combine_counts {

    tag { library.baseName }

    publishDir path: "${params.outputDir}/counts/${library.baseName}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file(counts) from countedFiles.collect()
    file(library) from combineLibraryFile

    output:
    file("counts_mageck.txt") into combinedMageckFile

    script:
    """
    combine_counts.R ${library} ${counts} > counts_mageck.txt
    """
}

fastqcSplitFiles = params.add_unknown_to_fastqc ? fastqcSplitFilesUnknown : fastqcSplitFiles

process fastqc {

    tag { id }

    input:
    set val(lane), val(id), file(fastq) from fastqcSplitFiles

    output:
    file "*_fastqc.{zip,html}" into fastqcResults

    script:
    """
    fastqc -t ${task.cpus} -q ${fastq}
    """
}

process multiqc {

    tag { 'all' }

    publishDir path: "${params.outputDir}",
               mode: 'copy',
               overwrite: 'true'

    input:
    file (fastqc: 'fastqc/*') from fastqcResults.collect()
    file (align: 'align/*') from alignResults.collect()
    file (featurecounts: 'featureCounts/*') from featureCountsResults.collect()

    output:
    file "*multiqc_report.html" into multiqc_report

    script:
    """
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    multiqc -f -x *.run .
    """
}

workflow.onComplete {
	println ( workflow.success ? "COMPLETED!" : "FAILED" )
}
