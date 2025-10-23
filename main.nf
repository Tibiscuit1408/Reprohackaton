nextflow.enable.dsl=2

process verify_docker {
    """
    echo "Packages installed :"
    fastq-dump --version
    bowtie --version
    samtools --version
    featureCounts -v
    R --version
    """
}

process fetch_data {
    tag {accession}
    publishDir "data", mode: 'copy'

    input:
    val accession

    output:
    path "data/*"

    script:
    """
    mkdir -p data
    fasterq-dump --progress --split-files --outdir data ${accession}
    """
}

workflow {
    verify_docker() 

    accessions = [
        "SRR10379721",
        "SRR10379722",
        "SRR10379723",
        "SRR10379724",
        "SRR10379725",
        "SRR10379726"
    ]

    Channel
        .from(accessions)
        .set { accession_ch }

    fetch_data(accession_ch)
}
