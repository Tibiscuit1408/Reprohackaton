#!/usr/bin/env nextflow

// 1. Generate the date-hour-min string
def run_timestamp = new Date().format('yyyyMMdd_HHmm')

// 2. Define your output parameter using the timestamp
params.outdir = "../Results/TEST_${run_timestamp}"


// 3. Trimming 
process Trimming{
    publishDir path: "${params.outdir}", mode: 'copy'
    input:
    path reads

    output:
    path "${reads.simpleName}_trimmed.fastq.gz"

    script:
    """
    cutadapt -m 25 -a TCAAGCAGAAGACGGCATACGAGATTACAAGGTGACTGGAGTTCAGACG \
        -o ${reads.simpleName}_trimmed.fastq.gz ${reads}
    """
}

process Mapping_refgenome{
    
    publishDir path: "${params.outdir}", mode: 'copy'
    input:
    path fasta

    output:
    path "Staphylococus.*"

    script:
    """
    bowtie-build ${fasta} Staphylococus
    """
}

process Mapping_replicate{
    
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path trimmed
    path index_files

    output:
    path "${trimmed.simpleName}_mapped.sam"

    script:
    """
    bowtie -x Staphylococus ${trimmed} -S ${trimmed.simpleName}_mapped.sam
    """
}

process Count{
    
    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path sam
    path gtf

    output:
    path "${sam.simpleName}_gene_counts.txt"

    script:
    """
    featureCounts -a ${gtf} -o ${sam.simpleName}_gene_counts.txt -g gene_id -t CDS -s 1 ${sam}
    """
}

workflow {
    reads_ch = Channel.fromPath('Test_head.fastq')
    fasta_ch = Channel.fromPath('reference_genome_Staphylococcus_aureus.fasta')
    gtf_ch   = Channel.fromPath('reference_genome_S_aureus.gtf')

    trimmed_ch = Trimming(reads_ch)
    index_ch   = Mapping_refgenome(fasta_ch)

    sam_ch = Mapping_replicate(trimmed_ch, index_ch)
    Count(sam_ch, gtf_ch)
}
