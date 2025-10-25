nextflow.enable.dsl=2

process fetch_data {
    container 'nf-sratoolkit'

    tag {accession}
    publishDir "data", mode: 'copy'

    input:
    val accession

    output:
    path "${accession}_1.fastq"

    script:
    """
    fasterq-dump --progress --outdir . ${accession}
    """
}

workflow {

    accession_ch = Channel
        .fromPath('accessions.txt')
        .splitText()
        .map { it.trim() }
        .filter { it }

    fetch_data(accession_ch)
}
