docker build -t reprohackaton:latest .
docker run -it --rm \
  -v "$(pwd):$(pwd)" \
  -w "$(pwd)"\
  reprohackaton:latest \
  bash -c "nextflow run main.nf"
