#!/usr/bin/env Rscript

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
    stop("Usage: Graphs.R <counts_file>")
}

counts_file <- args[1]

counts <- read.table(counts_file, header = TRUE, sep = "\t", comment.char = "#")

# Rename columns
colnames(counts)[colnames(counts) == "SRR10379726_1_trimmed_mapped.sam"] <- "Ctrl3"
colnames(counts)[colnames(counts) == "SRR10379723_1_trimmed_mapped.sam"] <- "IP3"
colnames(counts)[colnames(counts) == "SRR10379722_1_trimmed_mapped.sam"] <- "IP2"
colnames(counts)[colnames(counts) == "SRR10379724_1_trimmed_mapped.sam"] <- "Ctrl1"
colnames(counts)[colnames(counts) == "SRR10379725_1_trimmed_mapped.sam"] <- "Ctrl2"
colnames(counts)[colnames(counts) == "SRR10379721_1_trimmed_mapped.sam"] <- "IP1"

rownames(counts) <- counts[,1]
counts <- counts[,-1]

counts <- counts[, sapply(counts, is.numeric)]

counts_filtered <- counts[apply(counts, 1, var, na.rm=TRUE) != 0, ]
counts_filtered <- na.omit(counts_filtered)

sample_data <- counts[, c(4:ncol(counts_filtered))]

colData <- data.frame(
  condition = ifelse(grepl("^Ctrl", colnames(sample_data)), "control", "persister"),
  row.names = colnames(sample_data)
)

dds <- DESeqDataSetFromMatrix(
    countData = sample_data,
    colData = colData,
    design = ~ condition
)

dds <- DESeq(dds)

### VST transform
vsd <- vst(dds, blind = FALSE)

res <- results(dds)
resOrdered <- res[order(res$padj), ]
DEGs <- subset(resOrdered, padj < 0.1 & abs(log2FoldChange) > 1)

vst_mat <- assay(vsd)[rownames(vsd) %in% rownames(DEGs), ]
vst_mat_scaled <- t(scale(t(vst_mat)))

pdf("heatmap_DEGs.pdf", width = 8, height = 10)

annotation_col <- data.frame(Condition = factor(colData$condition))
rownames(annotation_col) <- rownames(colData)

my_colors <- colorRampPalette(c("deepskyblue4", "azure", "darkorange2"))(255)
ann_colors <- list(Condition = c(control = "black", persister = "grey"))

pheatmap(
    vst_mat_scaled,
    annotation_col = annotation_col,
    show_rownames = FALSE,
    show_colnames = FALSE,
    color = my_colors,
    annotation_colors = ann_colors,
    clustering_distance_rows = "correlation",
    clustering_distance_cols = "correlation",
    clustering_method = "ward.D2",
    main = "Heatmap of DEGs"
)

dev.off()

pdf("volcano_plot.pdf", width = 8, height = 6)

log_fc <- res$log2FoldChange
pval <- res$padj

valid_idx <- !is.na(log_fc) & !is.na(pval)
log_fc <- log_fc[valid_idx]
pval <- pval[valid_idx]

seuil <- 0.1
logfc_cutoff <- 1

plot(log_fc, -log10(pval),
     pch = 16,
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted p-value",
     main = "Volcano Plot",
     col = "lightgray")

abline(h = -log10(seuil), v = c(-logfc_cutoff, logfc_cutoff),
       lwd = 2, col = "orange")

sig_idx <- which(pval < seuil & abs(log_fc) > logfc_cutoff)
points(log_fc[sig_idx], -log10(pval[sig_idx]), pch = 16, col="red")

grid()

dev.off()
