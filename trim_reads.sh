
#!/bin/bash

mkdir trimmed_reads

read1=raw_reads/SRR031715_1.fastq.gz
read2=raw_reads/SRR031715_2.fastq.gz

OutputForwardPaired=trimmed_reads/SRR031715_1.trimmed.paired.fastq.gz
OutputForwardUnpaired=trimmed_reads/SRR031715_1.trimmed.unpaired.fastq.gz

OutputReversePaired=trimmed_reads/SRR031715_2.trimmed.paired.fastq.gz
OutputReverseUnpaired=trimmed_reads/SRR031715_2.trimmed.unpaired.fastq.gz

threads=8

java -jar Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads $threads $read1 $read2 \
$OutputForwardPaired $OutputForwardUnpaired \
$OutputReversePaired $OutputReverseUnpaired \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:true HEADCROP:3 TRAILING:10 MINLEN:25

