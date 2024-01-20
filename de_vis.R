#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

library(tidyverse) # V1.3.2
library(DESeq2) # V1.34.0
library(ggrepel) # V0.9.2
library(pheatmap) # V1.0.12
library(ggvenn) # V0.1.9
library(data.table) # V1.14.6

set.seed(1255342)

# Visualisation - volcano plot and heatmap

# 1. DGF
dds <- readRDS(file = "DGF_dds.rds")
res.volcano <- results(dds, alpha = 0.05)
res.volcano <- na.omit(as.data.frame(res.volcano))
sigGenes <- res.volcano %>% 
  filter(padj<0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  pull("gene")
res.volcano <- res.volcano %>% 
  mutate(sig = ifelse(padj < 0.05, "padj<0.05", "padj>0.05"))
lab <- res.volcano$padj < 0.05
res.volcano$lab <- row.names(res.volcano)
res.volcano$lab[!lab] <- NA
ggplot(data = res.volcano, aes(x=log2FoldChange, y=-log10(pvalue), colour=sig, label=lab)) +
  geom_point(size=2) +
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50')+
  guides(colour=guide_legend(title = NULL)) +
  labs(x= "Log2 Fold Change", y="-log10 (p value)")+
  geom_text_repel(size = 3, box.padding = 0.3)
res.heatmap <- as.data.frame(results(dds, alpha = 0.05))
res.heatmap <- rownames_to_column(res.heatmap, var = "gene")
sigGenes <- res.heatmap %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) 
sigGenes <- sigGenes %>% 
  pull("gene")
plotDat <- vst(dds, blind = FALSE)[sigGenes,] %>% 
  assay()
anno_info <- as.data.frame(colData(dds)[,c("Sample", "DGF_bin")])
anno_info <- anno_info[order(anno_info$DGF_bin),]
plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
anno_info$Sample <- NULL
plotDat2 <- t(scale(t(plotDat)))
trim <- function(x) {
  if (x > threshold){
    return(threshold)
  } else if (x < -threshold){
    return(-threshold)
  } else {
    return(x)
  }
} 
threshold <- 4
plotDat2 <-  apply(plotDat2, c(1,2), trim)
pheatmap(mat = plotDat, 
         scale="row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_info,
         gaps_col = 190,
         gaps_row = 9)
pheatmap(mat = plotDat2, 
         scale="none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_info,
         gaps_col = 190,
         gaps_row = 9)

# 2. EGFR3
dds <- readRDS(file = "EGFR3_dds.rds")
res.volcano <- results(dds, alpha = 0.05)
res.volcano <- na.omit(as.data.frame(res.volcano))
sigGenes <- res.volcano %>% 
  filter(padj<0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  pull("gene")
res.volcano <- res.volcano %>% 
  mutate(sig = ifelse(padj < 0.05, "padj<0.05", "padj>0.05"))
lab <- res.volcano$padj < 0.05
res.volcano$lab <- row.names(res.volcano)
res.volcano$lab[!lab] <- NA
ggplot(data = res.volcano, aes(x=log2FoldChange, y=-log10(pvalue), colour=sig, label=lab)) +
  geom_point(size=2) +
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50')+
  guides(colour=guide_legend(title = NULL)) +
  labs(x= "Log2 Fold Change", y="-log10 (p value)")+
  geom_text_repel(size = 3, box.padding = 0.3)
res.heatmap <- as.data.frame(results(dds, alpha = 0.05))
res.heatmap <- rownames_to_column(res.heatmap, var = "gene")
sigGenes <- res.heatmap %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) 
sigGenes <- sigGenes[c(1:20, (nrow(sigGenes)-19):nrow(sigGenes)),]
sigGenes <- sigGenes %>% 
  pull("gene")
plotDat <- vst(dds, blind = FALSE)[sigGenes,] %>% 
  assay()
anno_info <- as.data.frame(colData(dds)[,c("Sample", "EGFR3_cont")])
anno_info <- anno_info[order(anno_info$EGFR3_cont),]
plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
anno_info$Sample <- NULL
plotDat2 <- t(scale(t(plotDat)))
trim <- function(x) {
  if (x > threshold){
    return(threshold)
  } else if (x < -threshold){
    return(-threshold)
  } else {
    return(x)
  }
} 
threshold <- 3
plotDat2 <-  apply(plotDat2, c(1,2), trim)
pheatmap(mat = plotDat, 
         scale="row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_info)
pheatmap(mat = plotDat2, 
         scale="none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_info)

# 3. EGFR12
dds <- readRDS(file = "EGFR12_dds.rds")
res.volcano <- results(dds, alpha = 0.05)
res.volcano <- na.omit(as.data.frame(res.volcano))
sigGenes <- res.volcano %>% 
  filter(padj<0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  pull("gene")
res.volcano <- res.volcano %>% 
  mutate(sig = ifelse(padj < 0.05, "padj<0.05", "padj>0.05"))
lab <- res.volcano$padj < 0.05
res.volcano$lab <- row.names(res.volcano)
res.volcano$lab[!lab] <- NA
ggplot(data = res.volcano, aes(x=log2FoldChange, y=-log10(pvalue), colour=sig, label=lab)) +
  geom_point(size=2) +
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50')+
  guides(colour=guide_legend(title = NULL)) +
  labs(x= "Log2 Fold Change", y="-log10 (p value)")+
  geom_text_repel(size = 3, box.padding = 0.3)
res.heatmap <- as.data.frame(results(dds, alpha = 0.05))
res.heatmap <- rownames_to_column(res.heatmap, var = "gene")
sigGenes <- res.heatmap %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) 
sigGenes <- sigGenes[c(1:20, (nrow(sigGenes)-19):nrow(sigGenes)),]
sigGenes <- sigGenes %>% 
  pull("gene")
plotDat <- vst(dds, blind = FALSE)[sigGenes,] %>% 
  assay()
anno_info <- as.data.frame(colData(dds)[,c("Sample", "EGFR12")])
anno_info <- anno_info[order(anno_info$EGFR12),]
plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
anno_info$Sample <- NULL
summary(colnames(plotDat)==row.names(anno_info))
plotDat2 <- t(scale(t(plotDat)))
trim <- function(x) {
  if (x > threshold){
    return(threshold)
  } else if (x < -threshold){
    return(-threshold)
  } else {
    return(x)
  }
} 
threshold <- 4
plotDat2 <-  apply(plotDat2, c(1,2), trim)
pdf(file="de_analysis/plots/heatmap_EGFR12.pdf", useDingbats = FALSE)
pheatmap(mat = plotDat, 
         scale="row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_info,
         gaps_row = 20)
pheatmap(mat = plotDat2, 
         scale="none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_info,
         gaps_row = 20)