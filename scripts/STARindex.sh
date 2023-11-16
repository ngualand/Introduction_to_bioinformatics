#!/bin/bash

STAR --runMode genomeGenerate --runThreadN 1 --genomeDir "YOUR DIRECTORY FOR STORING GENOME" --genomeFastaFiles "GENOME FASTA FILE" --sjdbGTFfile "GTF ANNOTATION FILE" --sjdbOverhang "SPLICE SITES OVERHANG" --genomeSAindexNbases 11
