#!/usr/bin/env nextflow
process Featurecounts {

    publishDir "Data/featurecounts", mode: 'copy'

    input:
    path mapped_file

    output:
    path "${mapped_file.simpleName}_counts.txt"

    script:
    """
    featureCounts -a Data/annotation.gtf \
                  -o ${mapped_file.simpleName}_counts.txt \
                  $mapped_file
    """
}
