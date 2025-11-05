
# Projet Reprohackathon 

## Step 1 : Build the docker images 
```bash 
cd Dockerfile
docker build -t nf-bowtie -f Dockerfile.bowtie .
docker build -t nf-cutadapt -f Dockerfile.cutadapt .
docker build -t nf-samtools -f Dockerfile.samtools .
docker build -t nf-subread -f Dockerfile.subread .
docker build -t nf-sratoolkit -f Dockerfile.sratoolkit .
docker build -t gffread-ncbi -f Dockerfile.download . 
docker build -t nf-r-deseq -f Dockerfile.r-deseq .
```

## Step 2 : Launch Nextflow
```bash
cd ../Nextflow
nextflow run main.nf -with-docker -c main.config
```




