# Base Image
FROM python:3.8-slim-buster

# Metadata
LABEL base.image="python:3.8-slim-buster"
LABEL version="1.0"
LABEL software="CINTHIA"
LABEL software.version="202001"
LABEL description="an open source software tool to predict alpha-helical transmembrane topology"
LABEL website="https://github.com/BolognaBiocomp/cinthia"
LABEL documentation="https://github.com/BolognaBiocomp/cinthia"
LABEL license="GNU GENERAL PUBLIC LICENSE Version 3"
LABEL tags="Proteomics"
LABEL maintainer="Castrense Savojardo <castrense.savojardo2@unibo.it>"

WORKDIR /usr/src/cinthia

COPY requirements.txt .

RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt && \
    apt-get -y update && \
    apt-get install -y ncbi-blast+ wget && \
    useradd -m cinthia

WORKDIR /seqdb/

RUN wget https://share.biocomp.unibo.it/biocomp/sp2021_01/uniprot_sprot.fasta.gz && \
    gunzip uniprot_sprot.fasta.gz && \
    makeblastdb -in uniprot_sprot.fasta -dbtype prot

WORKDIR /usr/src/cinthia

USER cinthia

COPY . .

WORKDIR /data/

ENV CINTHIA_ROOT='/usr/src/cinthia' PATH=/usr/src/cinthia:$PATH

ENTRYPOINT ["/usr/src/cinthia/cinthia.py"]
