#!/bin/bash
cp -r /localdisk/data/BPSM/ICA1/fastq fastqseq
mkdir fastqc_output
cd fastqseq
gunzip *.gz
fastqc *.fq
unzip '*.zip'
mv *.fq ~/demo/fastqc_output
