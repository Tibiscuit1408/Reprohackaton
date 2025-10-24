docker build -t nf-nextflow -f Dockerfile.nextflow .
docker build -t nf-bowtie -f Dockerfile.bowtie .
docker build -t nf-cutadapt -f Dockerfile.cutadapt .
docker build -t nf-samtools -f Dockerfile.samtools .
docker build -t nf-subread -f Dockerfile.subread .
docker build -t nf-sratoolkit -f Dockerfile.sratoolkit .
nextflow run main.nf -with-docker
