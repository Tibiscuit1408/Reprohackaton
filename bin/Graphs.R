### Load required libraries

#Stat analysis (plots // article)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)

#Further analysis (for the report)
library(ggplot2)
library(reshape2)
library(dplyr)

############################
### Statistical Analysis ###
############################


### Import raw counts (reproduced data)
counts <- read.table("all_samples_gene_counts.txt",header = TRUE,sep = "\t",comment.char = "#")

#Rename the columns of the data
colnames(counts)[colnames(counts) == "SRR10379726_1_trimmed_mapped.sam"] <- "Ctrl3"
colnames(counts)[colnames(counts) == "SRR10379723_1_trimmed_mapped.sam"] <- "IP3"
colnames(counts)[colnames(counts) == "SRR10379722_1_trimmed_mapped.sam"] <- "IP2"
colnames(counts)[colnames(counts) == "SRR10379724_1_trimmed_mapped.sam"] <- "Ctrl1"
colnames(counts)[colnames(counts) == "SRR10379725_1_trimmed_mapped.sam"] <- "Ctrl2"
colnames(counts)[colnames(counts) == "SRR10379721_1_trimmed_mapped.sam"] <- "IP1"

#remove "gene-" of the gene_id so that it matches original data
counts$Geneid <- sub("^gene-", "", counts$Geneid)

#Set gene_id as the name of the rows
rownames(counts) <- counts[,1]

counts <- counts[,-1]

# Remove possible non-numeric columns
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

### Plot the heatmap ###
annotation_col <- data.frame(Condition = factor(colData$condition)) #aggregation of the controls and persisters (names of the )
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
         main = "Heatmap")


### Vulcano plot ###
# Volcano plot variables (from the DESeq analysis)
log_fc <- res$log2FoldChange       # log fold change
pval <- res$padj                   # adjusted p-value

# Optionally, remove NAs
valid_idx <- !is.na(log_fc) & !is.na(pval)
log_fc <- log_fc[valid_idx]
pval <- pval[valid_idx]

# Significance thresholds
seuil <- 0.1           # adjusted p-value threshold
logfc_cutoff <- 1      # log2 fold change threshold

# Basic volcano plot
plot(log_fc, -log10(pval),
     pch = 16,
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted p-value",
     main = "Volcano Plot",
     col = "lightgray")

# Add effect annotation
mtext("Effect: condition")

# Add significance lines
abline(h = -log10(seuil), v = c(-logfc_cutoff, logfc_cutoff),
       lwd = 2, col = "orange")

# Highlight significant DEGs
sig_idx <- which(pval < seuil & abs(log_fc) > logfc_cutoff)
points(log_fc[sig_idx], -log10(pval[sig_idx]), pch = 16)
# Optional grid
grid()

############################################################
### Difference between original data and reproduced data ###
############################################################

original <- read.table("GSE139659.tsv", sep='\t', header=TRUE)
## Rows and columns reorganisation

#set gene_id as name of the rows
rownames(original) <- original[,2]

# Remove possible non-numeric columns
original <- original[, sapply(original, is.numeric)]

original_filtered <- original[apply(original, 1, var, na.rm=TRUE) != 0, ]
original_filtered <- na.omit(original_filtered)

#plot heatmap

#select columns of interest
data <- original_filtered[,c(2:7)]

#reorganisation of rows based on ASC geneid
data <- data[order(rownames(data)), ]

#rename the columns
colnames(data)[colnames(data) == "ctrl4"] <- "Ctrl1"
colnames(data)[colnames(data) == "ctrl5"] <- "Ctrl2"
colnames(data)[colnames(data) == "ctrl6"] <- "Ctrl3"

#difference of length
#difference of genes selected and not selected

#plot of the regression -> amongst genes that are common in both datasets
## Global overview
df_merged <- merge(data, sample_data, by = "row.names")

# Identify columns by suffix
original_cols <- grep("\\.x$", names(df_merged), value = TRUE)
reproduced_cols <- gsub("\\.x$", ".y", original_cols)  # replace .x with .y to find matching names

#Combining columns original/repruced based on accurate column name
df_long <- do.call(rbind, lapply(seq_along(original_cols), function(i) {
  data.frame(
    Condition = gsub("\\.x$", "", original_cols[i]),
    Original = df_merged[[original_cols[i]]],
    Reproduced = df_merged[[reproduced_cols[i]]]
  )
}))

#Calculating the global MSE
df_long$SquaredError <- (df_long$Original - df_long$Reproduced)^2
MSE <- mean(df_long$SquaredError, na.rm = TRUE)
RMSE <- sqrt(MSE)

###normalised RMSE
nRMSE <- RMSE / mean(df_long$Original, na.rm = TRUE)
#nRMSE

nRMSE_per_condition <- df_long %>%
  group_by(Condition) %>%
  summarise(
    MSE = mean((Original - Reproduced)^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    Mean_Original = mean(Original, na.rm = TRUE),
    nRMSE = RMSE / Mean_Original
  )

#nRMSE_per_condition

### Statistics (how statistically different each pair is from each other)
p_values <- df_long %>%
  group_by(Condition) %>%
  summarise(
    p_value = t.test(Original, Reproduced, paired = TRUE)$p.value
  )

#p_values


### Regression (plot)
ggplot(df_long, aes(x = Original, y = Reproduced, color = Condition)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_text(
    data = nRMSE_per_condition,
    aes(x = Inf, y = Inf, label = paste("nRMSE =", round(nRMSE*100, 3),"%")),
    inherit.aes = FALSE,
    hjust = 2, vjust = 1.1,
    size = 3
  ) +
  geom_text(
    data = p_values,
    aes(x = Inf, y = Inf, label = paste("pvalue =",round(p_value,3))),
    inherit.aes = FALSE,
    hjust = 2.5, vjust = 3,
    size = 3
  ) +
  facet_wrap(~ Condition) +
  labs(
    title = "Regression Between Original and Reproduced Data",
    subtitle = "Green = perfect reproducibility (y = x)",
    x = "Original Data",
    y = "Reproduced Data"
  ) +
  theme_minimal()
