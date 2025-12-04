#!/usr/bin/env Rscript

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)
library(clusterProfiler)
library(KEGGREST)

#######################################################
### Analyse Statistique des données Reproduites ###
#######################################################



### Importer les donénes reproduites
counts_file <- args[1]

counts <- read.table(counts_file, header = TRUE, sep = "\t", comment.char = "#")

### Pré-traitements
#Renommer les colonnes
colnames(counts)[colnames(counts) == "SRR10379726_1_trimmed_mapped.bam"] <- "Ctrl3"
colnames(counts)[colnames(counts) == "SRR10379723_1_trimmed_mapped.bam"] <- "IP3"
colnames(counts)[colnames(counts) == "SRR10379722_1_trimmed_mapped.bam"] <- "IP2"
colnames(counts)[colnames(counts) == "SRR10379724_1_trimmed_mapped.bam"] <- "Ctrl1"
colnames(counts)[colnames(counts) == "SRR10379725_1_trimmed_mapped.bam"] <- "Ctrl2"
colnames(counts)[colnames(counts) == "SRR10379721_1_trimmed_mapped.bam"] <- "IP1"

rownames(counts) <- counts[,1]

counts <- counts[,-1]

# Remove possible non-numeric columns
counts <- counts[, sapply(counts, is.numeric)]

counts_filtered <- na.omit(counts)



# Sélection des 6 dernières colonnes
reproduced <- counts_filtered[,c(4:ncol(counts_filtered))]

### Define conditions automatically
colData <- data.frame(
  condition = ifelse(grepl("^Ctrl", colnames(reproduced)), "control", "persister"),
  row.names = colnames(reproduced)
)

colData$condition <- factor(colData$condition)

### Create DESeq2 dataset and run DE analysis
dds <- DESeqDataSetFromMatrix(countData = reproduced,
                              colData = colData,
                              design = ~ condition)

### Piepline DESeq2
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

## Extraction des gènes différentiellement exprimés
res <- results(dds)
resOrdered <- res[order(res$padj), ]
DEGs <- subset(resOrdered, padj < 0.1 & abs(log2FoldChange) > 1)
vst_mat <- assay(vsd)[rownames(vsd) %in% rownames(DEGs), ]
vst_mat_scaled <- t(scale(t(vst_mat)))

### Heatmap
annotation_col <- data.frame(Condition = factor(colData$condition))
rownames(annotation_col) <- rownames(colData)
my_colors <- colorRampPalette(c("deepskyblue4", "azure", "darkorange2"))(255)
ann_colors <- list(Condition = c(control = "black", persister = "grey"))
pdf("heatmap_vst_DEGs.pdf", width = 6, height = 6)
pheatmap(vst_mat_scaled,
         annotation_col = annotation_col,
         show_rownames = FALSE,
         show_colnames = FALSE,
         color = my_colors,
         annotation_colors = ann_colors,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         main = "")
dev.off()

### Vulcano plot

#Création d'un data frame supplémentaire ayant 1 colonne caractérisant l'état de régularisation des gènes
pdf("volcano_plot.pdf", width = 6, height = 6)
res_2 <- res
res_2$diffexpressed <- "NO"
  # if log2Foldchange > 1 and pvalue < 0.1, set as "UP"
  res_2$diffexpressed[res_2$log2FoldChange > 1 & res_2$padj < 0.1] <- "UP"
  # if log2Foldchange < -1 and pvalue < 0.1, set as "DOWN"
  res_2$diffexpressed[res_2$log2FoldChange < -1 & res_2$padj < 0.1] <- "DOWN"
  head(res_2[order(res_2$padj) & res_2$diffexpressed == 'DOWN', ])

res_2$gene_symbol <- rownames(res_2)
# Création d'une colonne "delabel" qui contiendra le nom des 15 gènes les plus différentiellement exprimés (NA sinon)
res_2$delabel <- ifelse(res_2$gene_symbol %in% head(res_2[order(res_2$padj), "gene_symbol"], 15), res_2$gene_symbol, NA)

#récupérer DEGs 
reg_repro <- data.frame(gene_symbol = res_2[!is.na(res_2$delabel), ]$gene_symbol,diffexpressed= res_2[!is.na(res_2$delabel), ]$diffexpressed)

ggplot(data = res_2, aes(x = log2FoldChange, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "white", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.1), col = "white", linetype = 'dashed') + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("#00AFBB", "grey", "#bb0c00"),  
                     labels = c("Sous-exprimé", "Non significatif", "Sur-exprimé")) + 
  coord_cartesian(ylim = c(0, 100), xlim = c(-10, 10)) + 
  labs(color = 'Niveau de régulation',  
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + 
  geom_text_repel(max.overlaps = Inf) 
dev.off()
############################################################################
### Différence entre les données de l'article et les données reproduites ###
############################################################################
original <- read.csv("GSE139659_IPvsctrl.csv",sep=";", header=TRUE)

#Nommer les lignes selon le nom des gènes
rownames(original) <- original[,2]

#Selectionner les colonnes d'intérêt
article <- original[,c(6:11)]
#Nommer les lignes selon nom des gènes
article <- article[order(rownames(article)), ]

#Renommer les colonnes
colnames(article)[colnames(article) == "ctrl4"] <- "Ctrl1"
colnames(article)[colnames(article) == "ctrl5"] <- "Ctrl2"
colnames(article)[colnames(article) == "ctrl6"] <- "Ctrl3"


#Analyse des gènes dans les données Articles et Reproduced

variables_analysis <- data.frame(
  Variable = union(rownames(article), rownames(reproduced)),
  In_Article = union(rownames(article), rownames(reproduced)) %in% rownames(article),
  In_Reproduced = union(rownames(article), rownames(reproduced)) %in% rownames(reproduced)
)
#Liste des gènes uniquement présents dans les données de l'article

genes_only_in_article <- variables_analysis %>%
  filter(In_Article == TRUE & In_Reproduced == FALSE) %>%
  pull(Variable)

#Listes des gènes uniquement présents dans les données reproduites
genes_only_in_reproduced <- variables_analysis %>%
  filter(In_Article == FALSE & In_Reproduced == TRUE) %>%
  pull(Variable)

#Regression (corrélation entre les deux datasets : Article et Reproduced)
df_merged <- merge(article, reproduced, by = "row.names")

# retirer toutes les lignes contenant des valeurs manquantes
df_merged <- na.omit(df_merged)



#Nommer les lignes selon les gènes
rownames(df_merged) <- df_merged$Row.names
df_merged$Row.names <- NULL

original_cols <- grep("\\.x$", names(df_merged), value = TRUE)
reproduced_cols <- gsub("\\.x$", ".y", original_cols)  # replace .x with .y to find matching names

df_long <- do.call(rbind, lapply(seq_along(original_cols), function(i) {
  data.frame(
    Gene = rownames(df_merged),                    
    Condition = gsub("\\.x$", "", original_cols[i]),
    Article = df_merged[[original_cols[i]]],
    Reproduced = df_merged[[reproduced_cols[i]]],
    stringsAsFactors = FALSE
  )
}))


#Calcul de la MSE globale et RMSE
df_long$SquaredError <- (df_long$Article - df_long$Reproduced)^2
MSE <- mean(df_long$SquaredError, na.rm = TRUE)
RMSE <- sqrt(MSE)

###Noramlisation de la RMSE
nRMSE <- RMSE / mean(df_long$Article, na.rm = TRUE)
#nRMSE

nRMSE_per_condition <- df_long %>%
  group_by(Condition) %>%
  summarise(
    MSE = mean((Article - Reproduced)^2, na.rm = TRUE),
    RMSE = sqrt(MSE),
    Mean_Article = mean(Article, na.rm = TRUE),
    nRMSE = RMSE / Mean_Article
  )

#nRMSE_per_condition

### Calcul des p_values

p_values <- df_long %>%
  group_by(Condition) %>%
  summarise(
    p_value = wilcox.test(Article, Reproduced, paired = TRUE)$p.value
  )
#p_values

pdf("regression_article_vs_reproduced.pdf", width = 7, height = 5)
### Regression (plot)
ggplot(df_long, aes(x = Article, y = Reproduced, color = Condition)) +
  geom_point(size = 2, alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, color = "darkgreen", linetype = "dashed") +
  facet_wrap(~ Condition) +
  labs(
    x = "Article",
    y = "Reproduction"
  ) +
  theme_minimal()
dev.off()

### graph Bland-Altman MA (Reproductibilité)

### Analyse Bland-Altman
## Préparation des données
#Récupération des log2FoldChange des données publiées et reproduites
prep_Degs_article <- original[,c("log2FoldChange","padj")]
prep_DEGs_article_2 <- subset(prep_Degs_article,padj < 0.1 & abs(log2FoldChange) > 1)
DEGs_article <- data.frame(Gene= rownames(prep_DEGs_article_2), log2FoldChange = prep_DEGs_article_2$log2FoldChange)
#Retirer la colonne Gene et la convertir en rownames
rownames(DEGs_article) <- DEGs_article$Gene
DEGs_article$Gene <- NULL
DEGs_reproduced <- as.data.frame(DEGs)["log2FoldChange"] #conversion d'un objet  “DESeqResults” en data frame
#length(rownames(original_filtered))
#length(rownames(prep_DEGs_article_2))
#length(rownames(DEGs_article))
#length(rownames(DEGs_reproduced))

DEGs_merged <- merge(DEGs_article, DEGs_reproduced, by = "row.names")
rownames(DEGs_merged) <- DEGs_merged$Row.names
DEGs_merged$Row.names <- NULL
degs_article_cols <- grep("\\.x$", names(DEGs_merged), value = TRUE)
degs_reproduced_cols <- gsub("\\.x$", ".y",degs_article_cols)


#Analyse des gènes dans les données Articles et Reproduced

DEGs_analysis <- data.frame(
  Variable = union(rownames(DEGs_article ), rownames(DEGs_reproduced)),
  In_Article = union(rownames(DEGs_article ), rownames(DEGs_reproduced)) %in% rownames(DEGs_article),
  In_Reproduced = union(rownames(DEGs_article ), rownames(DEGs_reproduced)) %in% rownames(DEGs_reproduced)
)
#Liste des gènes uniquement présents dans les données de l'article

DEGs_only_in_article <- DEGs_analysis %>%
  filter(In_Article == TRUE & In_Reproduced == FALSE) %>%
  pull(Variable)

#Listes des gènes uniquement présents dans les données reproduites
DEGs_only_in_reproduced <- DEGs_analysis %>%
  filter(In_Article == FALSE & In_Reproduced == TRUE) %>%
  pull(Variable)

df_long_2 <- do.call(rbind, lapply(seq_along(degs_article_cols), function(i) {
  data.frame(
    Gene = rownames(DEGs_merged),
    Condition = gsub("\\.x$", "", degs_article_cols[i]),
    Article = DEGs_merged[[degs_article_cols[i]]],
    Reproduced = DEGs_merged[[degs_reproduced_cols[i]]],
    stringsAsFactors = FALSE
  )
}))

df_ba_2 <- data.frame(
  Mean = (df_long_2$Article + df_long_2$Reproduced) / 2,
  Diff = df_long_2$Reproduced - df_long_2$Article,
  Condition = df_long_2$Condition
)

df_ba <- data.frame(
  Mean = (df_long_2$Article + df_long_2$Reproduced) / 2,
  Diff = df_long_2$Reproduced - df_long_2$Article,
  Condition = df_long_2$Condition
)
## Calcul des différences entre les logfold2change des données de l'Article et Reproduites
mean_diff <- mean(df_ba$Diff, na.rm = TRUE)
sd_diff   <- sd(df_ba$Diff, na.rm = TRUE)
LI <- mean_diff - 1.96 * sd_diff
LS <- mean_diff + 1.96 * sd_diff
df_ba$outlier_BA <- df_ba$Diff < LI | df_ba$Diff > LS
genes_significatifs_BA <- df_ba[df_ba$outlier_BA, ]

pdf("bland_altman_plot.pdf", width = 7, height = 5)
##Graphique Bland-Altman
df_ba$Gene <- rownames(DEGs_merged)
ggplot(df_ba, aes(x = Mean, y = Diff)) +
  geom_point(aes(color = outlier_BA), alpha = 0.6) +
  scale_color_manual(values = c("grey60", "red")) +
  geom_text_repel(
    data = df_ba[df_ba$outlier_BA == TRUE, ],
    aes(label = Gene),
    size = 3
  ) +
  geom_hline(yintercept = 0, color = "darkgreen", linetype = "dashed") +
  geom_hline(yintercept = mean(df_ba$Diff), color = "red")+
  geom_hline(yintercept = LI, color = "blue", linetype = "dotted") +
  geom_hline(yintercept = LS, color = "blue", linetype = "dotted") +
  theme_minimal()+
  labs(
    x = "Expression moyenne (Article + Reproduction)/2",
    y = "Différence (Reproduction − Article)",
    color = "Significativité"
  )
dev.off()

################################################
### Analyse Statistique des données Publiées ###
################################################


#Création d'un data frame supplémentaire ayant 1 colonne caractérisant l'état de régularisation des gènes
res_3 <- data.frame(log2FoldChange = original$log2FoldChange,
                    padj = original$padj,
                    gene_symbol = rownames(original))
res_3$diffexpressed <- "NO"
# if log2Foldchange > 1 and pvalue < 0.1, set as "UP"
res_3$diffexpressed[res_3$log2FoldChange > 1 & res_3$padj < 0.1] <- "UP"
# if log2Foldchange < -1 and pvalue < 0.1, set as "DOWN"
res_3$diffexpressed[res_3$log2FoldChange < -1 & res_3$padj < 0.1] <- "DOWN"
head(res_3[order(res_3$padj) & res_3$diffexpressed == 'UP', ])

# Création d'une colonne "delabel" qui contiendra le nom des 15 gènes les plus différentiellement exprimés (NA sinon)
res_3$delabel <- ifelse(res_3$gene_symbol %in% head(res_3[order(res_3$padj), "gene_symbol"], 15), res_3$gene_symbol, NA)

#Récupérer les lignes correspondant aux DEGs de l'article avec labels up/down
reg_article <- data.frame(gene_symbol = res_3[!is.na(res_3$delabel), ]$gene_symbol,diffexpressed= res_3[!is.na(res_3$delabel), ]$diffexpressed)


############################################################
####      MA-plot of all the RNASeq genes              #####
############################################################
# Set pdf document for output
pdf(file = "MA_plot.pdf", width = 5, height = 5) 

# Set colors for plot
col_black <- rgb(0, 0, 0, alpha = 0.3)     # transparent black for general points
col_red   <- rgb(1, 0, 0, alpha = 0.3)     # transparent red for significant points

# Identify extreme points (|log2foldchange| > 4) (& no NAs)
extreme         <- !is.na(res$log2FoldChange) & abs(res$log2FoldChange) > 4
extremely_big   <- !is.na(res$log2FoldChange) & res$log2FoldChange > 4
extremely_small <- !is.na(res$log2FoldChange) & res$log2FoldChange < -4

# Significant points
sig <- !is.na(res$padj) & res$padj < 0.05

# Base plot with only NON-extreme points
with(res[!extreme&!sig, ], plot(log10(baseMean), log2FoldChange,
                           pch = 21, cex = 0.3, col = col_black, bg = col_black,
                           xlab = "log10 mean of normalized counts",
                           ylab = "log2 fold change",
                           ylim = c(-4,4)))

abline(h = 0, col = 'black', lty = 4)

# Significant NON-extreme points
points(log10(res$baseMean)[sig & !extreme],
       res$log2FoldChange[sig & !extreme],
       col = col_red, bg = col_red, pch = 21, cex = 0.3)

# Significant extreme points
## Points with a log2foldchange below -4 are plotted at -4 as downward-pointing triangles 
y_extremely_small <- rep(-4, sum(sig & extremely_small))
points(log10(res$baseMean)[sig & extremely_small],
       y_extremely_small,
       col = col_red, bg = col_red, pch = 25, cex = 0.3)

## Points with a log2foldchange above 4 are plotted at 4 as upward-pointing triangles 
y_extremely_big   <- rep(4, sum(sig & extremely_big))
points(log10(res$baseMean)[sig & extremely_big],
       y_extremely_big,
       col = col_red, bg =col_red, pch = 24, cex = 0.3)

legend(x = -1, y=-3 ,
       legend = c("Significant", "Non-significant"), 
       col = c('red', 'grey'),          
       pch = c(16,16),
       cex=0.6)                

# Ouput
title("MA-plot of complete RNA-seq dataset")
dev.off()

############################################################
#####    RETRIEVING TRANSLATION RELATED GENES          #####
############################################################

# 1. Import the genes from kegg https://www.kegg.jp/kegg-bin/get_htext#C180 
kegg_path_id <- download_KEGG('sao')

## Import the translation related genes
id_ribosome <- kegg_path_id$KEGGPATHID2NAME$from[kegg_path_id$KEGGPATHID2NAME$to == 'Ribosome']
id_t_RNA <- kegg_path_id$KEGGPATHID2NAME$from[kegg_path_id$KEGGPATHID2NAME$to == 'Aminoacyl-tRNA biosynthesis']
gene_translation_path <- kegg_path_id$KEGGPATHID2EXTID$to[kegg_path_id$KEGGPATHID2EXTID$from %in% c(id_ribosome,id_t_RNA)]                                                     

## Add one gene from Ribosome biogenesis in eukaryotes that does not appear
gene_translation_path <- c(gene_translation_path, "SAOUHSC_01203")

## Import the genes coding for tRNA synthetases 
tRNA_synth_genes <- c("SAOUHSC_00509", "SAOUHSC_01722", "SAOUHSC_01737", 
                      "SAOUHSC_01471", "SAOUHSC_01666", "SAOUHSC_01788",
                      "SAOUHSC_00009", "SAOUHSC_00511", "SAOUHSC_00461",
                      "SAOUHSC_01183", "SAOUHSC_01767", "SAOUHSC_01875",
                      "SAOUHSC_01159", "SAOUHSC_00493", "SAOUHSC_00611",
                      "SAOUHSC_01240", "SAOUHSC_01738", "SAOUHSC_01092",
                      "SAOUHSC_01093", "SAOUHSC_01839", "SAOUHSC_00933")

# 2. Import the genes from kegg's orthology cf. https://www.kegg.jp/brite/ko03012
ko_nums <- c("ko:K02518","ko:K02519","ko:K02520","ko:K02357","ko:K02358",
             "ko:K02355","ko:K02356","ko:K03833","ko:K04568","ko:K19810",
             "ko:K09906","ko:K02835","ko:K02836","ko:K02837","ko:K02838",
             "ko:K02839","ko:K09890","ko:K15034","ko:K02493","ko:K01056",
             "ko:K04794")
kegg_res = lapply(ko_nums, keggLink, target="genes") #may take some time
kegg_ortho_id = data.frame(value = unname(unlist(kegg_res)))

gene_translation_ortho <- c()
for (gene_id in kegg_ortho_id$value) {
  if (substring(gene_id,1,4) == 'sao:') {
    gene_translation_ortho <- c(gene_translation_ortho, substring(gene_id,5))
  }
}

# 3. Add genes from NCBI (not found in kegg)
## gene coding for the ribosome recycling factor
frr_gene <- 'SAOUHSC_01236'

## gene coding for the peptidyl-tRNA hydrolase
pth_gene <- 'SAOUHSC_00475'

## 1rst gene coding for the translation initiation factor 
infA_gene <- 'SAOUHSC_02489'

## 2nd gene coding for the translation initiation factor 
infB_gene <- 'SAOUHSC_01246'

## 3rd gene coding for the translation initiation factor 
infC_gene <- 'SAOUHSC_01786'

## 4 genes coding for the translation elongation factor
tsf_genes <- c('SAOUHSC_00529', 'SAOUHSC_00530', 'SAOUHSC_01234', 'SAOUHSC_01625')

# sum(c(frr_gene, pth_gene, infA_gene, infB_gene, infC_gene, tsf_genes) %in% gene_translation_ortho) 
# these genes are already in gene_translation_ortho

# 4. Final vector with all the retrieved translation related genes 
gene_translation <- c(gene_translation_path, gene_translation_ortho)

#length(gene_translation) # 170
#sum(gene_translation %in% rownames(res)) #170
#sum(is.na(res[gene_translation,"baseMean"])) #0 all genes have a calculated baseMean
#sum(is.na(res[gene_translation,"log2FoldChange"])) #1 gene without log2FoldChange
#sum(is.na(res[gene_translation,"padj"])) #2 genes with no adjusted p-value 
#sum(!is.na(res[gene_translation,"padj"]) &  res[gene_translation,"padj"]< 0.05) #56 significant genes
# sum(!is.na(res[gene_translation,"padj"]) &  res[gene_translation,"padj"]> 0.05) #112
############################################################
####      MA-plot of genes related to translation      #####
############################################################

# Set pdf document for output
pdf(file = "MA_plot_translation.pdf", width = 5, height = 5) 

# Plot non-significant points from translation_genes
with(res[rownames(res) %in% gene_translation & !sig, ], 
     plot(log2(baseMean), log2FoldChange,
          pch = 21, cex = 0.6, col = col_black, bg = col_black,
          xlab = "log2 mean of normalized counts",
          ylab = "log2 fold change",
          xlim = c(0,20), ylim = c(-6,5), new=T))

# Plot significant points from translation_genes
with(res[rownames(res) %in% gene_translation & sig, ], 
     points(log2(baseMean), log2FoldChange,
            pch = 21, cex = 0.6, col = 'red', bg = 'red'))

# Plot significant genes coding for t-RNA synthetases
with(res[rownames(res) %in% tRNA_synth_genes, ], 
     points(log2(baseMean), log2FoldChange,
            pch = 1, cex = 0.6, col = 'black'))

## frr
with(res[frr_gene, ], 
     arrows(x0=log2(baseMean), y0=log2FoldChange,
            x1=10, y1=2, col='black', lwd=1, angle=0))
text(x=9.7, y=2.3, labels="frr", col='black')

## pth
with(res[pth_gene, ], 
     arrows(x0=log2(baseMean), y0=log2FoldChange,
            x1=7, y1=-3, col='black', lwd=1, angle=0))
text(x=6.7, y=-3.3, labels="pth", col='black')

## infA
with(res[infA_gene, ], 
     arrows(x0=log2(baseMean), y0=log2FoldChange,
            x1=12, y1=2.7, col='black', lwd=1, angle=0))
text(x=12, y=3, labels="infA", col='black')

## infB
with(res[infB_gene, ], 
     arrows(x0=log2(baseMean), y0=log2FoldChange,
            x1=17, y1=0.5, col='black', lwd=1, angle=0))
text(x=17.9, y=0.5, labels="infB", col='black')

## infC
with(res[infC_gene, ], 
     arrows(x0=log2(baseMean), y0=log2FoldChange,
            x1=17, y1=1, col='black', lwd=1, angle=0))
text(x=17.9, y=1, labels="infC", col='black')

#tsf gene
with(res[tsf_genes[3], ], 
     arrows(x0=log2(baseMean), y0=log2FoldChange,
            x1=17, y1=3.5, col='black', lwd=1, angle=0))
text(x=17, y=3.8, labels="tsf", col='black')

legend(x = 0, y=-4 ,
       legend = c("Significant", "Non-significant","AA-tRNA synthetases"), 
       col = c('red', 'grey', 'black'),          
       pch = c(16,16,1),
       cex=0.6)                

# Output
abline(h = 0, col = 'black', lty = 4)
title("MA-plot of genes related to translation")

dev.off()
