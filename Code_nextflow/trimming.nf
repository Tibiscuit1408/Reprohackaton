#!/usr/bin/env nextflow

process Trimming {

    publishDir "Data", mode: 'copy'

    input:
    path input_file

    output:
    path "${input_file.simpleName}_trimmed.fastq"

    script:
    """
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
             -m 25 \
             -o ${input_file.simpleName}_trimmed.fastq \
             $input_file
    """
}
