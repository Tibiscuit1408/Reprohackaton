library(clusterProfiler)
library(KEGGREST)

############################################################
####      MA-plot of all the RNASeq genes              #####
############################################################
# Set pdf document for output
pdf(file = "/bin/MA_plot.pdf", width = 5, height = 5) 

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
kegg_tr_genes <- read.table("C:/Users/iness/Documents/Post BAC/AgroParisTech/3eme_annee/AMI2B/Reprohackathon/genes_related_to_translation.txt",sep = "\t")
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
pdf(file = "/bin/MA_plot_translation.pdf", width = 5, height = 5) 

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