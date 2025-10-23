FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV PATH="/usr/local/bin:/usr/local/bin/sratoolkit/sratoolkit.3.2.1-ubuntu64/bin:${PATH}"

# Update and install of all dependencies
RUN apt-get update --fix-missing
RUN apt-get install -y \
    wget curl git unzip build-essential \
    r-base \
    openjdk-17-jdk \
    bowtie \
    samtools \
    subread \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/local/bin

# Nextflow
RUN wget -qO nextflow https://get.nextflow.io && \
    chmod +x nextflow

# SRA Toolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.current-ubuntu64.tar.gz && \
    mkdir -p sratoolkit && \
    mv sratoolkit.* sratoolkit && \
    ln -s /usr/local/bin/sratoolkit/bin/* /usr/local/bin/ && \
    rm -f sratoolkit.current-ubuntu64.tar.gz

CMD ["bash"]
