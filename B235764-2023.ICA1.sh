#!/bin/bash
mkdir ICA-1
cd ICA-1
cp -r /localdisk/data/BPSM/ICA1/fastq fastqseq
mkdir fastqc_output
mkdir Countsdata
cd fastqseq
gunzip *.gz
fastqc *.fq
unzip '*.zip'
mv ~/ICA-1/fastqseq/*.fq ~/ICA-1/fastqc_output

# This first bit copies the raw data files to a directory fastqseq. This is then followed by an extraction of the .gz files and performing a quality check on them using fastqc.After unzipping the outputs of the fastqc program, the resultant .fq files are moved to a seperate folder 'fastqc_output', thus finishing the first task

cd ..
cp -r /localdisk/data/BPSM/ICA1/Tcongo_genome ~/ICA-1/Alignment
cd Alignment
gunzip *.gz
cd ..
cd ~/ICA-1/fastqc_output
cp -r ~/ICA-1/fastqc_output/*.fq ~/ICA-1/Alignment
cd ..
cd Alignment

# The Tcongo genome file and all the .fq files are then copied on to a seperate folder to perform the alignment.The genome gz file is extracted. We now have all the required files ready for the alignment

bowtie2-build ~/ICA-1/Alignment/TriTrypDB-46_TcongolenseIL3000_2019_Genome.fasta IL_3000_index

# Using bowtie2 command to create an index for the Tcongolense fasta file for the alignent. Since the.fasta format is not in the bowtie index format, the bowtie2 code written above creates an index for the same.

index_prefix="IL_3000_index"
for file in ~/ICA-1/Alignment/*_1.fq; do
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

#This for/do section loops through all files with the .fq extention. It places the file name as variables and replaces _.fq with _2.fq when it comes across those files. The loop then proceeds to allign all files individually with the reference IL-3000 genome index generated using bowtie2  and provides the resultant output in sam format. This output file  is then read and indexed by samtools in a bam format for the alignment of each read pair with the reference IL-3000 genome index. This is the end and completion of Task 3

cd ..
cp -r /localdisk/data/BPSM/ICA1/TriTrypDB-46_TcongolenseIL3000_2019.bed ~/ICA-1/Countsdata
cd Countsdata
cp -r ~/ICA-1/Alignment/*_output.bam .
cp -r ~/ICA-1/Alignment/*_output.bam.bai .
for file in *.bam
        {
                bedtools coverage -counts -a TriTrypDB-46_TcongolenseIL3000_2019.bed -b ${file} > ${file}.txt
        }

# This bit copies the .bed file of TcongolenseIL3000 along with the output.bam and it's index files to the same directory. When all three file types are present, a simple for loop, loops through all the .bam files with the .bed file and creates respective output.txt files which provide the count data for each of the file. This is the completion of Task 4.

cp -r ~/ICA-1/Alignment/*.txt ~/ICA-1/Average

# We now copy all the .txt files with the counts data generated in the previous task into a new folder

cd ..
cd Average

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

# Manually referring to the Tco2.fqfiles from the first step we write code to group the 48 individual samples into 15 different groups by merging all the sample txt files in a particular group with one another

cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
cat CL1_0_Unin.txt | cut -f6 | awk '{ sum += $1 } END { if (NR > 0) print "Average:",sum/NR}'
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

# We now display the files and note down the average for each group. The above script cuts out the 6th collumn with the number of alignments and proceeds to create an average of the entire collumn using an awk script which adds the values of all integers in the sixth collumn and divides it by the number of lines (entries). We then note down the values for all 15 groups

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

# This final bit creates a .tsv file with 4 collumns detailing the Sample type, Time, Condition and the Average for the respective group. Although a loop could have been made to simplify this process, i have performed the same manually due to lesser number of samples and groups. I recognise that this method is not feasible with larger sets and amounts of data. It then displays the tab separated output file. This is the completion of Task 5.

# Additional comments
# I have tried to use simple linux commands and scripts to run the tasks. I have not made any changes to the program or tool parameters.

