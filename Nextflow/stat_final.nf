#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def run_timestamp = new Date().format('yyyyMMdd_HHmm')
params.outdir = "../Results/TEST_${run_timestamp}"

process Rstat {

    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path counts_file

    output:
    path "*.pdf"

    script:
    """
    Graphs_mercredi.R ${counts_file}
    """
}

workflow {

    counts_ch = Channel.fromPath("count_gene.txt")

    Rstat(counts_ch)
}
