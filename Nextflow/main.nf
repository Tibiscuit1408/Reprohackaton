#!/usr/bin/env nextflow
nextflow.enable.dsl=2
// 1. Generate a timestamp for unique output folder names
def run_timestamp = new Date().format('yyyyMMdd_HHmm')

// 2. Define your output directory parameter
params.outdir = "../Results/TEST_${run_timestamp}"

process fetch_data {
    publishDir "${params.outdir}/data", mode: 'copy'

    input:
    val accession

    output:
    path("${accession}*.fastq")

    script:
    """
    fasterq-dump --progress --split-files --outdir . ${accession}
    """
}
process get_genomic_data {

    publishDir "./Genomic_data", mode: 'symlink'

    output:
    path "reference_genome_Staphylococcus_aureus.fasta"
    path "reference_genome_S_aureus.gtf"

    script:
    """
    # Download the dataset
    datasets download genome accession GCF_000013425.1 --include genome,gff3 --filename ncbi_dataset.zip

    # Unzip it
    unzip -o ncbi_dataset.zip

    # Move FASTA and GFF files to current directory
    mv ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna reference_genome_Staphylococcus_aureus.fasta
    mv ncbi_dataset/data/GCF_000013425.1/genomic.gff reference_genome_S_aureus.gff

    # Convert GFF to GTF
    gffread reference_genome_S_aureus.gff -T -o reference_genome_S_aureus.gtf
    """
}

// 3. Trimming 
process Trimming {

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

// 4. Build Bowtie index
process Mapping_refgenome {

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

// 5. Map reads to reference
process Mapping_replicate {

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

// 6. Count reads per gene
process Count {

    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path sams
    path gtf

    output:
    path "all_samples_gene_counts.txt"

    script:
    """
    featureCounts -a ${gtf} -o all_samples_gene_counts.txt -g transcript_id -t CDS -s 1 ${sams.join(' ')}
    """
}


workflow {
    accession_ch = Channel
        .fromPath('accessions.txt')
        .splitText()
        .map { it.trim() }
        .filter { it }
    reads_ch = fetch_data(accession_ch)

    (fasta_ch, gtf_ch) = get_genomic_data()
    

    trimmed_ch = Trimming(reads_ch)
    index_ch   = Mapping_refgenome(fasta_ch)

    sam_ch = Mapping_replicate(trimmed_ch, index_ch)
    sam_list_ch = sam_ch.collect()
    Count(sam_list_ch, gtf_ch)
}
