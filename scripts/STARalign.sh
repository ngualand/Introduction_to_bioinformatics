#!/bin/bash

STAR --genomeDir "STARindex" --readFilesCommand gunzip -c --runThreadN 1 --readFilesIn "FASTQ1 and FASTQ2" --outFileNamePrefix "FILENAMEPREFIX" --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outSAMattributes All --quantMode GeneCounts
