library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyverse)
library(Factoshiny)


#######################################################
### Analyse Statistique sur les données Reproduites ###
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
         main = "")


### Vulcano plot
# Volcano plot variables (from the DESeq analysis)
log_fc <- res$log2FoldChange       # log fold change
pval <- res$padj                   # adjusted p-value

valid_idx <- !is.na(log_fc) & !is.na(pval)
log_fc <- log_fc[valid_idx]
pval <- pval[valid_idx]

# Seuils
seuil <- 0.1           # adjusted p-value threshold
logfc_cutoff <- 1      # log2 fold change threshold

plot(log_fc, -log10(pval),
     pch = 16,
     xlab = "log2 Fold Change",
     ylab = "-log10 Adjusted p-value",
     main = "",
     col = "lightgray")

mtext("Effect: condition")
abline(h = -log10(seuil), v = c(-logfc_cutoff, logfc_cutoff),
       lwd = 2, col = "orange")
sig_idx <- which(pval < seuil & abs(log_fc) > logfc_cutoff)
points(log_fc[sig_idx], -log10(pval[sig_idx]), pch = 16)
# Optional grid
grid()

############################################################################
### Différence entre les données de l'article et les données reproduites ###
############################################################################

original <- read.table("GSE139659.tsv", sep='\t', header=TRUE)

#Nommer les lignes selon le nom des gènes
rownames(original) <- original[,2]

# Convertir les données en valeurs numériques
original <- original[, sapply(original, is.numeric)]
# retirer toutes les lignes contenant des valeurs manquantes (30 lignes supprimées)
original_filtered <- na.omit(original)


#Selectionner les colonnes d'intérêt
article <- original_filtered[,c(2:7)]

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
    p_value = t.test(Article, Reproduced, paired = TRUE)$p.value
  )

#p_values


### Regression (plot)
ggplot(df_long, aes(x = Article, y = Reproduced, color = Condition)) +
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
    x = "Article Data",
    y = "Reproduced Data"
  ) +
  theme_minimal()


### graph Bland-Altman MA (Reproductibilité)

### Analyse Bland-Altman
## Préparation des données
#Récupération des log2FoldChange des données publiées et reproduites
prep_Degs_article <- original_filtered[,c("log2FoldChange","padj")]
prep_DEGs_article_2 <- subset(prep_Degs_article,padj < 0.1 & abs(log2FoldChange) > 1)
DEGs_article <- data.frame(Gene= rownames(prep_DEGs_article_2), log2FoldChange = prep_DEGs_article_2$log2FoldChange)
#Retirer la colonne Gene et la convertir en rownames
rownames(DEGs_article) <- DEGs_article$Gene
DEGs_article$Gene <- NULL
DEGs_reproduced <- as.data.frame(DEGs)["log2FoldChange"] #conversion d'un objet  “DESeqResults” en data frame
#length(rownames(original_filtered))
#1988
#length(rownames(prep_DEGs_article_2))
#827
#length(rownames(DEGs_article))
#827
#length(rownames(DEGs_reproduced))
#1174

DEGs_merged <- merge(DEGs_article, DEGs_reproduced, by = "row.names")
rownames(DEGs_merged) <- DEGs_merged$Row.names
DEGs_merged$Row.names <- NULL
degs_article_cols <- grep("\\.x$", names(DEGs_merged), value = TRUE)
degs_reproduced_cols <- gsub("\\.x$", ".y",degs_article_cols)

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



### ACP et HCPC

select1 <- ((original_filtered$padj < 0.1) & (abs(original_filtered$log2FoldChange)>=1))
class(select1)
sum(select1)

genePositive_article = original_filtered[select1,]
genePositive_article <- original_filtered[,c(2:7)]

## Données de l'Article
x = as.matrix(genePositive_article)
dta = data.frame(t(x))
res.PCA<-PCA(dta,graph=TRUE)
res.HCPC <- HCPC(res.PCA, graph=FALSE)
#plot.HCPC(res.HCPC,choice='tree',title='Arbre hiérarchique')
plot.HCPC(res.HCPC,choice='map',draw.tree=FALSE,title='')
#plot.HCPC(res.HCPC,choice='3D.map',ind.names=FALSE,centers.plot=FALSE,angle=60,title='Arbre hiérarchique sur le plan factoriel')


desc_cluster2 = res.HCPC$desc.var$quanti[[2]]
desc_cluster3 = res.HCPC$desc.var$quanti[[3]]

## Données reproduites
genePositive_reproduced = counts_filtered[select1,]
head(genePositive_reproduced)
genePositive_reproduced <- counts_filtered[,c(4:9)]

x_2 = as.matrix(genePositive_reproduced)
dta_2 = data.frame(t(x_2))
res.PCA_2<-PCA(dta_2,graph=TRUE)
res.HCPC_2 <- HCPC(res.PCA_2, graph=FALSE)
#plot.HCPC(res.HCPC_2,choice='tree',title='Arbre hiérarchique')
plot.HCPC(res.HCPC_2,choice='map',draw.tree=FALSE,title='')
#plot.HCPC(res.HCPC_2,choice='3D.map',ind.names=FALSE,centers.plot=FALSE,angle=60,title='Arbre hiérarchique sur le plan factoriel')
# Distribution in clusters
clusters_2 = res.HCPC_2$data.clust$clust
#length(dta$c.phenotype.lignee.)
#length(clusters)

# Association gènes clusters  (focus on cluster 2 : Perissters)
desc_cluster2_2 = res.HCPC_2$desc.var$quanti[[2]]
head(desc_cluster2_2)
negative_associations_2_2 = rownames(desc_cluster2_2)[desc_cluster2_2[,1]<0]
positive_associations_2_2 = rownames(desc_cluster2_2)[desc_cluster2_2[,1]>0]

# Association gènes clusters (focus on cluster 3 : Contrôles)
desc_cluster3_2 = res.HCPC_2$desc.var$quanti[[3]]
#head(desc_cluster3_2)
negative_associations_3_2 = rownames(desc_cluster3_2)[desc_cluster3_2[,1]<0]
positive_associations_3_2 = rownames(desc_cluster3_2)[desc_cluster3_2[,1]>0]
