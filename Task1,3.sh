#!/bin/bash
cp -r /localdisk/data/BPSM/ICA1/fastq fastqseq
mkdir fastqc_output
cd fastqseq
gunzip *.gz
fastqc *.fq
unzip '*.zip'
mv *.fq ~/demo/fastqc_output
cd ..
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome Alignment
cd Alignment
gunzip *.gz
cd ..
cd fastqc_output
cp -r *.fq ~/demo/Alignment
cd ..
cd Alignment
bowtie2-build ~/demo/Alignment/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta IL_3000_index
index_prefix="IL_3000_index"
for file in ~/demo/Alignment/*_1.fq; do
    if [ -e "$file" ]; then
        filename1=$(basename "$file")
        filename2="${filename1/_1.fq/_2.fq}"
        output_sam="${filename1/_1.fq/_output.sam}"
        output_bam="${filename1/_1.fq/_output.bam}"
        bowtie2 -x "$index_prefix" -1 "$file" -2 "${file/_1.fq/_2.fq}" -S "$output_sam"
        samtools view -bS "$output_sam" | samtools sort -o "$output_bam"
        samtools index "$output_bam"
    fi
done

