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
    path("reference_genome_Staphylococcus_aureus.fasta")
    path("reference_genome_S_aureus.gtf")

    script:
    """
    # Download the dataset
    datasets download genome accession GCF_000013425.1 --include genome,gtf --filename ncbi_dataset.zip

    # Unzip it
    unzip -o ncbi_dataset.zip

    # Move FASTA and GFF files to current directory
    mv ncbi_dataset/data/GCF_000013425.1/GCF_000013425.1_ASM1342v1_genomic.fna reference_genome_Staphylococcus_aureus.fasta
    mv ncbi_dataset/data/GCF_000013425.1/genomic.gtf reference_genome_S_aureus.gtf
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
process Index_creation {

    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path fasta

    output:
    path "Staphylococus.*.ebwt"

    script:
    """
    bowtie-build ${fasta} Staphylococus
    """
}

// 5. Mapping reads and generating sorted BAM
process Mapping_replicate {

    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path trimmed
    path index_files

    output:
    path "${trimmed.simpleName}_mapped.bam"

    script:
    """
    # Map reads with Bowtie, convert SAM to BAM, and sort in one go
    zcat ${trimmed} | \
    bowtie -S Staphylococus - | \
    samtools view -bS - | \
    samtools sort -o ${trimmed.simpleName}_mapped.bam
    """
}


// 6. Count reads per gene using sorted BAM files
process Count {

    publishDir path: "${params.outdir}", mode: 'copy'

    input:
    path bams
    path gtf

    output:
    path "all_samples_gene_counts.txt"

    script:
    """
    featureCounts -a ${gtf} -o all_samples_gene_counts.txt -g gene_id -t gene -s 1 ${bams.join(' ')}
    """
}
//7. Statistic with R 
process Rstat {

    publishDir path: "${params.outdir}/stats", mode: 'copy'

    input:
    path counts_file

    output:
    path *.pdf

    script:
    """
    Graphs_mercredi.R ${counts_file}
    """
}



workflow {
    accession_ch = Channel
        .fromPath('Nextflow/accessions.txt')
        .splitText()
        .map { it.trim() }
    reads_ch = fetch_data(accession_ch)

    (fasta_ch, gtf_ch) = get_genomic_data()
    

    trimmed_ch = Trimming(reads_ch)
    index_ch   = Index_creation(fasta_ch)

    bam_ch = Mapping_replicate(trimmed_ch, index_ch)
    bam_list_ch = bam_ch.collect()
    count_file = Count(bam_list_ch, gtf_ch)
    Rstat(counts_file)
}



