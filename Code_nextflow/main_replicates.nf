#!/usr/bin/env nextflow

include { Trimming } from './trimming.nf'
include { Mapping } from './mapping.nf'
include { Featurecounts} from './featurecounts.nf'

Channel
    .fromPath("Source/persister/*.fastq.gz")
    .set { fastq_files }



workflow{
    Trimming(fastq_files)
    Mapping(Trimming.out)
    Featurecounts(Mapping.out)
    

}
