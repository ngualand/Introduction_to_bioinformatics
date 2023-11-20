#Install Libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

#import libraries
library(DESeq2)

#Import table to convert Ensembl ID to genename (build it via BioMart for example)
genes.names = read.delim("/home/nicolo/EnsemblGenes_name.bed", header = F)

# Read sample information
sample_info <- read.delim("/home/nicolo/Downloads/filereport_read_run_PRJNA305863_tsv.txt", header = TRUE, stringsAsFactors = FALSE)
#Set rown
rownames(sample_info) = sample_info$run_accession
sample_info$study_accession = sample_info$run_accession


# Create DESeqDataSet
dds <- DESeqDataSetFromHTSeqCount(sampleTable = sample_info,
                                  directory = "/home/nicolo/Documents/Uniud/Corso_Bioinfo/Files/htseqcounts",
                                  design = ~ Group)


# Perform DESeq analysis
dds <- DESeq(dds)

# Extract differential expression results
res <- results(dds, contrast = c("Group", "H2O2", "CTRL"))
res.genes = as.data.frame(res)
res.genes = merge(res.genes, genes.names, by.x = 0 , by.y = "V4", all.x = T)

# Save results to a CSV file
write.csv(res, file = "deseq_results.csv")

# Create MA plot
plotMA(res, main = "DESeq2 MA Plot")

#Transform data for PCA and heatmap
vsd = varianceStabilizingTransformation(dds, blind = T)

#Plot expression of Genes of interest
plotCounts(dds, gene=which.min(res$padj), intgroup="Group")

# Create a heatmap of the top differentially expressed genes
library(pheatmap)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$Group
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(4, "Blues")) )(4)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


#plot PCA
DESeq2::plotPCA(vsd, intgroup = c("Group"))


#Heatmpa of DE genes 
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:500]
df <- as.data.frame(colData(dds)[,c("Group")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=T, annotation_col= sample_info[,c("run_accession", "Group")])

#Find upregulated and downregulated genes 
up.genes = res.genes[res.genes$log2FoldChange > 0 & res.genes$padj < 0.1 & !is.na(res.genes$padj),]
down.genes = res.genes[res.genes$log2FoldChange < 0 & res.genes$padj < 0.1 & !is.na(res.genes$padj),]


#Enrichment analysis
library(gprofiler2)
#Upregulated
gene_list <- unlist(lapply(up.genes$Row.names, function(x) unlist(strsplit(x, "\\."))[1]))
background = unlist(lapply(res.genes$Row.names, function(x) unlist(strsplit(x, "\\."))[1]))

# Specify the organism and database you want to use
organism = "hsapiens"  # Human


# Perform the enrichment analysis
gostres = gost(gene_list,
               organism = organism,
               ordered_query = F,
               significant = TRUE, 
               exclude_iea = TRUE, 
               user_threshold = 0.1, 
               domain_scope = "custom",
               custom_bg = background, 
               evcodes = F,
               correction_method = "fdr", 
               sources = c("GO:BP", "REAC", "KEGG"))

xxx = gostres$result

p = gostplot(gostres, capped = TRUE, interactive = F)
publish_gostplot(p, highlight_terms = c("REAC:R-HSA-381042"), 
                 width = NA, height = NA, filename = NULL )

