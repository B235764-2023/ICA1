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
cd ..
cp -r /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed ~/demo/Countsdata
cd demo/Countsdata
cp -r ~/demo/Alignment/*_output.bam .
cp -r ~/demo/Alignment/*_output.bam.bai .
for file in *.bam
        {
                bedtools coverage -counts -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${file} > ${file}.txt
        }
cp -r *.txt ~/demo/Average
cat Tco-290_output.bam.txt Tco-907_output.bam.txt Tco-106_output.bam.txt Tco-580_output.bam.txt > CL2_0_Unin.txt
cat Tco-230_output.bam.txt Tco-486_output.bam.txt Tco-674_output.bam.txt > CL1_0_Unin.txt
cat Tco-719_output.bam.txt Tco-467_output.bam.txt Tco-392_output.bam.txt Tco-964_output.bam.txt > CL2_24_Unin.txt
cat Tco-503_output.bam.txt Tco-702_output.bam.txt Tco-960_output.bam.txt > CL2_24_in.txt
cat Tco-549_output.bam.txt Tco-483_output.bam.txt Tco-580_output.bam.txt > CL2_48_Unin.txt
cat Tco-622_output.bam.txt Tco-596_output.bam.txt Tco-582_output.bam.txt > CL2_48_in.txt
cat Tco-21_output.bam.txt Tco-859_output.bam.txt Tco-86_output.bam.txt > CL1_24_Unin.txt
cat Tco-935_output.bam.txt Tco-613_output.bam.txt Tco-851_output.bam.txt > CL1_24_in.txt
cat Tco-486_output.bam.txt Tco-999_output.bam.txt Tco-522_output.bam.txt Tco-362_output.bam.txt > CL1_48_Unin.txt
cat Tco-76_output.bam.txt Tco-28_output.bam.txt Tco-892_output.bam.txt Tco-397_output.bam.txt > CL1_48_in.txt
cat Tco-754_output.bam.txt Tco-229_output.bam.txt Tco-878_output.bam.txt > WT_24_Unin.txt
cat Tco-752_output.bam.txt Tco-633_output.bam.txt Tco-480_output.bam.txt Tco-427_output.bam.txt > WT_24_in.txt
cat Tco-17_output.bam.txt Tco-122_output.bam.txt Tco-159_output.bam.txt > WT_48_Unin.txt
cat Tco-444_output.bam.txt Tco-398_output.bam.txt Tco-480_output.bam.txt Tco-757_output.bam.txt > WT_48_in.txt
cat Tco-160_output.bam.txt Tco-949_output.bam.txt Tco-182_output.bam.txt > WT_0_Unin.txt
cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:", sum / NR }'
cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum / NR }'
cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR }'
cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL1_24_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL1_24_in.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL1_48_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL1_48_in.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL2_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL2_24_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL2_24_in.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL2_48_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL2_48_in.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat WT_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat WT_24_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat WT_24_in.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat WT_48_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat WT_48_in.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'

echo -e "Sampletype\tTime\tCondition\tAverage" > Averagedata.tsv
echo -e "Clone1\t0\tUninduced\t9.75304" >> Averagedata.tsv
echo -e "Clone1\t24\tUninduced\t10.1984" >> Averagedata.tsv
echo -e "Clone1\t24\tinduced\t10.6561" >> Averagedata.tsv
echo -e "Clone1\t48\tUninduced\t10.5592" >> Averagedata.tsv
echo -e "Clone1\t48\tinduced\t10.8104" >> Averagedata.tsv
echo -e "Clone2\t0\tUninduced\t9.74023" >> Averagedata.tsv
echo -e "Clone2\t24\tUninduced\t10.6347" >> Averagedata.tsv
echo -e "Clone2\t24\tinduced\t10.095" >> Averagedata.tsv
echo -e "Clone2\t48\tUninduced\t9.82383" >> Averagedata.tsv
echo -e "Clone2\t48\tinduced\t10.452" >> Averagedata.tsv
echo -e "WT\t0\tUninduced\t10.9588" >> Averagedata.tsv
echo -e "WT\t24\tUninduced\t11.2263" >> Averagedata.tsv
echo -e "WT\t24\tinduced\t10.9998" >> Averagedata.tsv
echo -e "WT\t48\tUninduced\t7.68739" >> Averagedata.tsv
echo -e "WT\t48\tinduced\t11.6171" >> Averagedata.tsv

cat Averagedata.tsv

