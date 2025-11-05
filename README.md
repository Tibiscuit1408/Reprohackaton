

##Step 1 - build the dockers:

cd Dockerfile
docker build -t nf-bowtie -f Dockerfile.bowtie .
docker build -t nf-cutadapt -f Dockerfile.cutadapt .
docker build -t nf-samtools -f Dockerfile.samtools .
docker build -t nf-subread -f Dockerfile.subread .
docker build -t nf-sratoolkit -f Dockerfile.sratoolkit .
docker build -t gffread-ncbi -f Dockerfile.download . 
#docker build -t nf-r-deseq -f Dockerfile.r-deseq .
#le dernier est long à installer et pas nécessaire pour l'instant (c'est seulement pour du R)

##Step 2:
#techniquement, le fetch data viendra avant
cd ../Nextflow
nextflow run main.nf -with-docker -c main.config
#nextflow run fetch_data.nf -with-docker -c fetch_data.config
#Fetchdata => ça ne marche pas, ça bug jsp pourquoi. 
#fastrq-dump est reconnu dans le docker 
#Rentrer dans le docker pour vérifier les commandes : 
#docker run -it nf-sratoolkit:latest
#Lancer ce nf là fait buguer mon ordi, il y a peut-être des histoires de cpu/mémoire à gérer
#la progression est censée être montrée aussi mais il y a un problème avec et 
#comme mon ordi plante à chaque fois, je sais pas quoi



