FROM rust:1.98.1

## Install cmake 
RUN apt update && apt install cmake -y

## Pull release 
RUN wget https://github.com/willhooper/bam-zero-mapq/archive/refs/tags/v0.1.tar.gz && tar -xvzf v0.1.tar.gz

## Install Riker
RUN cd bam-zero-mapq-0.1 && cargo build --release