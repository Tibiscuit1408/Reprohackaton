# Projet Reprohackathon
Ce projet a été réalisé afin d'être reproductible facilement. Les images des dockers ont été publiées pour faciliter le lancement du code nextflow. Voici l>

Toutes les images sont hébergées sur le dépôt Docker Hub suivant :
[https://hub.docker.com/repositories/tibiscuit14/reprohackaton2025](https://hub.docker.com/repositories/tibiscuit14/reprohackaton2025)

| Image utilisée dans Nextflow | Docker Hub URL |
|------------------------------|----------------|
| [`nf-bowtie`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=nf-bowtie) | `tibiscuit14/reprohackaton2025:nf-bowtie` |
| [`nf-cutadapt`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=nf-cutadapt) | `tibiscuit14/reprohackaton2025:nf-cutadapt` |
| [`nf-subread`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=nf-subread) | `tibiscuit14/reprohackaton2025:nf-subread` |
| [`nf-sratoolkit`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=nf-sratoolkit) | `tibiscuit14/reprohackaton2025:nf-sratoolkit` |
| [`gffread-ncbi`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=gffread-ncbi) | `tibiscuit14/reprohackaton2025:gffread-ncbi` |
| [`nf-r-deseq`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=nf-r-deseq) | `tibiscuit14/reprohackaton2025:nf-r-deseq` |
| [`nf-r-enrichmentbrowser`](https://hub.docker.com/r/tibiscuit14/reprohackaton2025/tags?name=nf-r-enrichmentbrowser) | `tibiscuit14/reprohackaton2025:nf-r-enrichmentbrowser` |

## OPTIONNAL - Etape 0 : Création des images dockers
Si vous souhaitez créer les images dockers sur votre ordinateur, voici les commandes nécessaires. Il faudra cependant changer les configurations du nextflo>
```bash
cd Dockerfile
docker build -t nf-bowtie -f Dockerfile.bowtie .
docker build -t nf-cutadapt -f Dockerfile.cutadapt .
docker build -t nf-subread -f Dockerfile.subread .
docker build -t nf-sratoolkit -f Dockerfile.sratoolkit .
docker build -t gffread-ncbi -f Dockerfile.download .
docker build -t nf-r-deseq -f Dockerfile.r-deseq .
```

## Etape 1 : Lancement de Nextflow
```bash
nextflow run ./Nextflow/main.nf -with-docker -c ./Nextflow/main.config
```
