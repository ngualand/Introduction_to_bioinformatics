#!/bin/bash

htseq-count -a 10 -s no -f bam -r pos --nonunique all  -m union --type exon -i gene_id BAMFILE GTFFILE > OUTPUT
