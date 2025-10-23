#!/usr/bin/env nextflow

process Mapping_ref_genome {

    publishDir "Data", mode: 'copy'

    input:
    path genome_fasta

    output:
    path "reference_genome_Staphylococcus_aureus_indexed.*"

    script:
    """
    bowtie-build $genome_fasta reference_genome_Staphylococcus_aureus_indexed
    """
}

params.genome = "Source/ref/reference_genome_Staphylococcus_aureus.fasta"

workflow {
    Mapping_ref_genome(params.genome)
}