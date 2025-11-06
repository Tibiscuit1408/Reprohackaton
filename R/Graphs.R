### Load required libraries
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

### Import raw counts
#counts <- read.table("GSE139659.csv", header=TRUE, sep=";")
counts <- read.table("all_samples_gene_counts.txt",header = TRUE,sep = "\t",comment.char = "#")

#Rename the columns of the data
colnames(counts)[colnames(counts) == "SRR10379726_1_trimmed_mapped.sam"] <- "Ctrl3"
colnames(counts)[colnames(counts) == "SRR10379723_1_trimmed_mapped.sam"] <- "IP3"
colnames(counts)[colnames(counts) == "SRR10379722_1_trimmed_mapped.sam"] <- "IP2"
colnames(counts)[colnames(counts) == "SRR10379724_1_trimmed_mapped.sam"] <- "Ctrl1"
colnames(counts)[colnames(counts) == "SRR10379725_1_trimmed_mapped.sam"] <- "Ctrl2"
colnames(counts)[colnames(counts) == "SRR10379721_1_trimmed_mapped.sam"] <- "IP1"

#Set gene_id as the name of the rows
rownames(counts) <- counts[,1]
counts <- counts[,-1]

# Remove possible non-numeric columns (safety)
counts <- counts[, sapply(counts, is.numeric)]

counts_filtered <- counts[apply(counts, 1, var, na.rm=TRUE) != 0, ]
counts_filtered <- na.omit(counts_filtered)



# Sélection des 6 dernières colonnes
sample_data <- counts[,c(4:ncol(counts_filtered))]

### Define conditions automatically
colData <- data.frame(
  condition = ifelse(grepl("^Ctrl", colnames(sample_data)), "control", "persister"),
  row.names = colnames(sample_data)
)

### Create DESeq2 dataset and run DE analysis
dds <- DESeqDataSetFromMatrix(countData = sample_data,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

### Apply variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

### Extraction of top DEGs
res <- results(dds)
resOrdered <- res[order(res$padj), ]
DEGs <- subset(resOrdered, padj < 0.1 & abs(log2FoldChange) > 1)

# Keep only genes in the VST data
vst_mat <- assay(vsd)[rownames(vsd) %in% rownames(DEGs), ]

### Center and scale per gene (optional but recommended)
vst_mat_scaled <- t(scale(t(vst_mat)))

### Plot the heatmap
annotation_col <- data.frame(Condition = factor(colData$condition)) #aggreagation of the controls and persisters (names of the )
rownames(annotation_col) <- rownames(colData)
my_colors <- colorRampPalette(c("deepskyblue4", "azure", "darkorange2"))(255)
ann_colors <- list(Condition = c(control = "black", persister = "grey"))
pheatmap(vst_mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = my_colors,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         main = NULL)


