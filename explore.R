#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

library(pheatmap) # V1.0.12
library(tidyverse) # V1.3.2
library(DESeq2) # V1.34.0
library(ggstatsplot) # V0.10.0
library(ggExtra) # V0.10.0

set.seed(1255342)

# key data exploration including clustering samples and PCA

# Import data
coldata <- read.csv("Coldata.csv")
row.names(coldata) <- coldata$Sample
countTable <- read.csv("CountTable.csv")
row.names(countTable) <- countTable$X
countTable$X <- NULL

# 1. Demonstrating importance of cortex/medulla composition in transcriptional profile
# Score based on ratio of ABCA13 (Med) : CRYAA (Cortex)
tmp <- as.data.frame(t(countTable[row.names(countTable) %in% c("ABCA13", "CRYAA"),]))
tmp$logratio <- log2((tmp$ABCA13+1)/(tmp$CRYAA+1))
ggplot(tmp, aes(x=row.names(tmp), y=(logratio)))+
  geom_point()+
  geom_hline(yintercept = 0)
coldata <- merge(coldata, tmp, by=0)
row.names(coldata) <- coldata$Row.names
coldata$Row.names <- NULL
coldata <- coldata[,!colnames(coldata)%in%c("ABCA13", "CRYAA")] # remove ABCA13 and CRYAA cols
coldata$CvMlogRatio <- coldata$logratio
coldata$logratio <- NULL
# Principal Component Analysis
data <- as.matrix(countTable)
data <- vst(data)
data <- data[order(rowVars(data), decreasing = TRUE),]
data <- data[c(1:500), ]
data <- t(data)
pca <- prcomp(data)
var_scores <- pca$sdev^2
var_scores <- var_scores/sum(var_scores)*100
barplot(var_scores, main = "Scree Plot", xlab = "Principal Component", ylab = "% Variation")
pc_scores <- data.frame(pca$x, coldata[row.names(coldata) %in% row.names(pca$x),])
ggplot(data = pc_scores, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = CvMlogRatio, shape = DSEX), size = 4) +
  geom_hline(yintercept = 0, color = 'grey70') +
  geom_vline(xintercept = 0, color = 'grey70') +
  xlab(paste("PC1 (", round(var_scores[1], 2), "%)", sep = "")) +
  ylab(paste("PC2 (", round(var_scores[2], 2), "%)", sep = "")) +
  ggtitle("PCA")
# PC1 = Cortex vs Medulla; PC2 = Sex
# Use scaled PC1 scores to account for tissue composition
# Calculate scaled PC1 scores
coldata <- merge(coldata, pc_scores[,c(1,2)], by=0)
row.names(coldata) <- coldata$Row.names
coldata$Row.names <- NULL
coldata$PC1Scaled <- scale(coldata$PC1) 


# 2. Cluster cortical and medullary samples for WGCNA analysis
# Import cortexmedulla genelist
genelist <- read.csv(file = "cortexmedulla_genelist.csv", header = TRUE, sep = ",")
sigGenes <- genelist$X
coldata$Sample_Set <- as.factor(coldata$Sample_Set)
model <- as.formula( ~ Sample_Set + DSEX + PC1Scaled)
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = model)
plotDat <- vst(dds, blind = FALSE)[sigGenes,] %>% 
  assay()
anno_info <- as.data.frame(colData(dds)[,c("PC1Scaled","Sample")])
plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
anno_info$PC1Scaled <- NULL
(out <- pheatmap(mat = plotDat,
                scale="row", 
                cluster_rows = FALSE, 
                cluster_cols = TRUE, 
                clustering_distance_cols="euclidean",
                show_rownames = TRUE,
                show_colnames = FALSE,
                annotation_col = anno_info))
colnames(countTable[,out$tree_col[["order"]]])
tmp <- as.data.frame(sort(cutree(out$tree_col, k=3)))
tmp <- rownames_to_column(tmp, var = "Sample")
tmp <- tmp[match(row.names(coldata), tmp$Sample),]
tmp <- mutate(tmp, Cluster = tmp[,2])
tmp <-  tmp[, c(1, 3)]
coldata <- merge(coldata, tmp, by="Sample")
coldata[,'Cluster']<-factor(coldata[,'Cluster'])
coldata <- droplevels(coldata)
ggplot(coldata, aes(x = Cluster, y = CvMlogRatio))+
  geom_jitter(aes(colour = Cluster))
# cluster 1 = medulla; cluster 3 = cortex; cluster 2 = mix
coldata <- coldata %>% 
  mutate(Tissue = ifelse(Cluster == 1,
                         "medulla",
                         ifelse(Cluster == 2,
                                "mixed",
                                "cortex")))
rownames(coldata) <- coldata$Sample


# 3. Evidence for Sample_Set batch effect
model <- as.formula(~ DSEX + PC1Scaled + Sample_Set)
reduced.model <- as.formula (~ DSEX + PC1Scaled)
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = model)
dds <- dds [rowSums(counts(dds)) >0 ,]
dds <- DESeq(dds, test = "LRT", reduced = reduced.model, betaPrior = FALSE)
res <- results(dds, alpha = 0.05)
sum(res$padj < 0.05, na.rm = TRUE)
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
  labs(x= "Log2 Fold Change", y="-log10 (p value)")

# save updated coldata and countTable files
write.csv(coldata, file = "Coldata.csv")
write.csv(countTable, file = "CountTable.csv")
