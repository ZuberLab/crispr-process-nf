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
                                        - sample_name:  name of demultiplexed sample followed by an integer number (1-4)
                                                        corresponding to the length of the stagger sequence used for this sample,
                                                        seperated by the keyword STAGGERLENGTH (e.g. STAGGERLENGTH1, STAGGERLENGTH2, ...)
                                        - barcode:      8 nucleotide long sequence of the stagger, followed by the sample barcode,
                                                        filled up to a length of 8 nucleotides with the subsequent sequence of the spacer

        --reverse_complement        Should reads be converted to reverse_complement
                                    before trimming (default: false)

        --barcode_random_length     Number of nucleotides in random barcode
                                    (default for sgRNA: 2)

        --barcode_demux_mismatches  Number of mismatches allowed during demultiplexing
                                    of barcode. (default: 2)

        --barcode_length            Number of nucleotides in sample barcode.
                                    (default: 4)

        --spacer_length             Number of nucleotides in spacer sequence between
                                    barcodes and sgRNA / shRNA sequence. (default: 20)

        --guide_length              Number of nucleotides in guide sequence. (default: 21)

        --padding_beginning         Nucleotides used for 5' padding if sgRNA / shRNA are of
                                    unequal length. Must be one of G, C, T, and A.
                                    (default: ACC)

        --padding_end               Nucleotides used for 3' padding if sgRNA / shRNA are of
                                    unequal length. Must be one of G, C, T, and A.
                                    (default: GGT)

     Profiles:
        standard                    local execution
        singularity                 local execution with singularity
        ii2                         SLURM execution with singularity on IMPIMBA2

     Docker:
     zuberlab/crispr-nf:latest

     Author:
     Jesse J. Lipp (jesse.lipp@imp.ac.at)
     Florian Andersch (florian.andersch@imp.ac.at)
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
log.info " input directory          : ${params.inputDir}"
log.info " output directory         : ${params.outputDir}"
log.info " library file             : ${params.library}"
log.info " barcode file             : ${params.barcodes}"
log.info " barcode random (nt)      : ${params.barcode_random_length}"
log.info " barcode demultiplex (nt) : ${params.barcode_length}"
log.info " spacer (nt)              : ${params.spacer_length}"
log.info " demultiplex mismatches   : ${params.barcode_demux_mismatches}"
log.info " 5' guide padding base    : ${params.padding_beginning}"
log.info " 3' guide padding base    : ${params.padding_end}"
log.info " reverse complement       : ${params.reverse_complement}"
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

senseFastqFiles = Channel.create()
antisenseFastqFiles = Channel.create()

fastqFilesFromBam
    .mix(rawFastqFiles)
    .choice( senseFastqFiles, antisenseFastqFiles ) { params.reverse_complement ? 1 : 0 }

process reverse_complement_fastq {

    tag { lane }

    input:
    set val(lane), file(fastq) from antisenseFastqFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into rcfastqFiles

    script:
    """
    mv ${fastq} input.fastq.gz
    zcat input.fastq.gz | fastx_reverse_complement \
        -z \
        -o ${lane}.fastq.gz
    """
}

senseFastqFiles
    .mix(rcfastqFiles)
    .set { fastqFiles }

process trim_random_barcode {

    tag { lane }

    input:
    set val(lane), file(fastq) from fastqFiles

    output:
    set val(lane), file("${lane}.fastq.gz") into randomBarcodeTrimmedFiles

    script:
    position = params.barcode_random_length
    """
    if [[ ${position} -gt 0 ]]
    then
      mv ${fastq} input.fastq.gz
      cutadapt -u ${position} -j ${task.cpus} input.fastq.gz -o ${lane}.fastq.gz
    else
      mv ${fastq} input.fastq.gz
      path_to_file=\$(readlink -f input.fastq.gz)
      cp \${path_to_file} ${lane}.fastq.gz
    fi
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
    .filter { it[1].baseName =~ /^(?!.*unknown).*$/ }
    .map { lane, file -> tuple(lane, file.name.replaceAll(/\.fq\.gz/, ''), file) }
    .into { flattenedSplitFiles; fastqcSplitFiles }

process trim_barcode_and_spacer {

    tag { id }

    publishDir path: "${params.outputDir}/fastq/${lane}",
               mode: 'copy',
               overwrite: 'true'

    input:
    set val(lane), val(id), file(fastq) from flattenedSplitFiles

    output:
    set val(lane), val(id), file("${id}.fastq.gz") into spacerTrimmedFiles

    script:
    barcode_spacer_length = params.spacer_length + params.barcode_length
    """
    str=${id}
    stagger_length="\$( awk -F 'STAGGERLENGTH' '{print \$2}' <<< \${str})"
    remove_beginning=\$(expr \${stagger_length} + ${barcode_spacer_length})

    cutadapt ${fastq} -j ${task.cpus} -u \${remove_beginning} -o ${id}_remove_beginning.fastq.gz
    cutadapt ${id}_remove_beginning.fastq.gz -j ${task.cpus} -l ${params.guide_length} -o ${id}.fastq.gz

    rm ${id}_remove_beginning.fastq.gz

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
    process_library.R ${library} ${params.padding_beginning} ${params.padding_end}
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
    """
    bowtie2 \
        --threads \$((${task.cpus})) \
        -x ${index}/index \
        -L ${params.guide_length} \
        --score-min 'C,0,-1' \
        -N 0 \
        --seed 42 \
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


fastqcSplitFiles
    .map { lane, id, file -> tuple(lane, file) }
    .groupTuple()
    .set { fastqcSplitFiles }

process fastqc {

    tag { lane }

    input:
    set val(lane), file(fastq) from fastqcSplitFiles

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
