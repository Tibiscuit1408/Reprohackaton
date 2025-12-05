#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def run_timestamp = new Date().format('yyyyMMdd_HHmm')
params.outdir = "../Results/TEST_${run_timestamp}"

process Rstat {

    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path counts_file
    path csv_or

    output:
    path "*.pdf"

    script:
    """
    Graphs_code.R ${counts_file} ${csv_or}
    """
}

workflow {

    counts_ch = Channel.fromPath("count_gene.txt")
    csv_or = Channel.fromPath("GSE139659_IPvsctrl.csv")

    Rstat(counts_ch,csv_or)
}
