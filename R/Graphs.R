library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggrepel)


#######################################################
### Analyse Statistique des données Reproduites ###
#######################################################


### Importer les donénes reproduites
counts <- read.table("count_gene.txt",header = TRUE,sep = "\t",comment.char = "#")

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


### Vulcano plot

#Création d'un data frame supplémentaire ayant 1 colonne caractérisant l'état de régularisation des gènes
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

