#!/usr/bin/env nextflow
process Mapping {

    publishDir "Data", mode: 'copy'

    input:
    path trimmed_file

    output:
    path "${input_file.simpleName}_mapped.sam"

    script:
    """
    bowtie -x reference_genome_Staphylococcus_aureus_indexed \
           -U $trimmed_file \
           -S ${trimmed_file.simpleName}_mapped.sam
    """
}
