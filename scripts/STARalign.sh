#!/bin/bash

#Map reads to the reference genome
STAR --genomeDir "STARindex" --readFilesCommand gunzip -c --runThreadN 1 --readFilesIn "FASTQ1 and FASTQ2" --outFileNamePrefix "FILENAMEPREFIX" --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattributes All --quantMode GeneCounts

#Create index for BAM files
samtools index BAMFILE
