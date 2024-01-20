## Analysis of RNASeq data from Kidney RNASeq of corex vs medulla. Samples from April 2016 batch.
## Analysis by J R Ferdinand comencing May 2016 using DESeq2 for DE.


setwd(".")
load(".Rdata")

library("DESeq2")
library("ggplot2")
library(VennDiagram)
library(RColorBrewer)
library(gplots)
library(pheatmap)


## Plot graphs from sequencing QC data provided by AROS
seqqc <- read.csv("sequencing_QC.csv", row.names = 1)
colnames(seqqc)<- c("qc30","reads","quality")  
  ## QC30
pdf(file = "QC30.pdf")
ggplot(data =seqqc, aes(x=row.names(seqqc), y=qc30)) +
  geom_bar(stat='identity', fill="deepskyblue1") +
  geom_hline(yintercept=0, color='black') +
  geom_hline(yintercept=80, color='black', linetype=3, size =2) +
  xlab ("Sample") +
  ylab ( expression("% bases" >= "QC30")  ) +
  ylim (0,100) +
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100))
dev.off()
  ##reads
seqqc$readsxmil <- seqqc$reads / 1000000

pdf(file = "reads.pdf")
ggplot(data =seqqc, aes(x=row.names(seqqc), y=readsxmil)) +
  geom_bar(stat='identity', fill="deepskyblue1") +
  geom_hline(yintercept=0, color='black') +
  xlab ("Sample") +
  ylab (expression("reads x 10"^"6")) +
  ylim (0,120)
dev.off()  

load("gene_counts_table.RData")
##Remove anotation apart from Gene name which should be in row.names
countTable <- counts_tab$counts
## Rename columns
colnames(countTable) <- c(1,10,2,3,4,5,6,7,8,9)
## Read in table of sample information located in CWD
coldata <- read.csv("RNAseqdesign.csv", row.names = 1)
## Re order columns
countTable <- countTable [, c(1,3,4,5,6,7,8,9,10,2)]
## Construct a DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = ~Person + Tissue)
## pre filter to remove rows with total count of less than 1
dds <- dds [ rowSums(counts(dds)) >1 ,]

## set medulla to be the refrence level
dds$Tissue <- relevel(dds$Tissue, ref="medulla")

## Differential expression analysis
dds <- DESeq(dds, betaPrior = T)

## Extract resualts of DE and print sumary to screen
DEres <- results(dds, alpha=0.05)
summary(DEres)
DiffExpression <- as.data.frame(DEres)
names(DiffExpression)[names(DiffExpression)=="log2FoldChange"] <- "log2FoldChange (cortex over medulla)"
write.csv( DiffExpression, file="DE_Cortex_vs_Medulla.csv", row.names = TRUE)

##QC plots

##MA plot
pdf(file="MA_plot.pdf")
plotMA(DEres, main="MA Plot", ylim=c(-6,6))
dev.off()
##Box plot of counts
## Add 0.1 to all counts
countTableqc <- lapply(countTable, function(x){x+0.1})
pdf(file="boxplot_of_raw_counts.pdf")
boxplot(countTableqc, las=2, ylab="Counts",  log="y", main="Box plot of raw counts")
dev.off()

## Transformation of data for visulisation 
## Regularised log transformation
rld <- rlog(dds, blind=FALSE)
## Variance stabilised transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
## Output normalized counts from VSD
write.csv(assay(vsd), file="VSD_counts.csv")


## PCA plot
pdf(file="PCA1_vs_PCA2_sample_groups.pdf")
plotPCA(rld, intgroup="Tissue")
dev.off()
pdf(file="PCA1_vs_PCA2_indervidules.pdf")
plotPCA(rld, intgroup="Person")
dev.off()

## PCA plot with prcomp
data <- assay(vsd)
data <- t(data) ## Need samples in rows and genes (or other variables) in columns
pca <- prcomp( data)
pcasum <- (summary(pca))$importance
scores <-data.frame(pca$x, coldata)
  ##PC1 vs 2
pdf(file="PCA_1v2_paired.pdf")
ggplot( data = scores, aes(x=PC1, y=PC2, color=Tissue, group= Person, shape=Sex)) +
  geom_point(size=8) +
  geom_line( colour="black", linetype = 3, size=1) +
  geom_hline(yintercept=0, color='grey70') +
  geom_vline(xintercept = 0, color='grey70') +
  xlab("PC1 (57.63% Variance)") +
  ylab("PC2 (12.24% Variance)")
dev.off()

##PC2 vs 3
pdf(file="PCA_2v3_paired.pdf")
ggplot( data = scores, aes(x=PC2, y=PC3, color=Tissue, group= Person, shape=Sex)) +
  geom_point(size=8) +
  geom_line( colour="black", linetype = 3, size=1) +
  geom_hline(yintercept=0, color='grey70') +
  geom_vline(xintercept = 0, color='grey70') +
  ylab("PC3 (8.465% Variance)") +
  xlab("PC2 (12.24% Variance)")
dev.off()


## Volcano plot of DE
names(DiffExpression)[names(DiffExpression)=="log2FoldChange (cortex over medulla)"] <- "log2FoldChange"


DiffExpression$sig <- as.factor(DiffExpression$padj < 0.05)
DiffExpressionNA <- na.omit(DiffExpression)
pdf(file="Volcano_Cortex_vs_Medulla.pdf")
ggplot(data = DiffExpressionNA, aes(x=log2FoldChange, y=-log10(pvalue), colour=sig)) +
  geom_point(size=1) +
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50')+
  guides(colour=guide_legend(title = NULL)) +
  scale_colour_manual(breaks = c("FALSE","TRUE"), values=c("grey33","red"), labels=c(expression("padj" >= "0.05"), "padj < 0.05"))
dev.off()

## Extreem points labled
DElab <- DiffExpressionNA
DElab$Gene <- row.names(DElab)
lab <-  -log10(DElab$pvalue)>25 &DElab$log2FoldChange>1 | -log10(DElab$pvalue)>38 &DElab$log2FoldChange<0 | -log10(DElab$pvalue)>20 &DElab$log2FoldChange<(-4)
DElab$Gene[!lab] <- NA

pdf(file="Volcano_Cortex_vs_Medulla_Labled.pdf")


ggplot(data = DElab, aes(x=log2FoldChange, y=-log10(pvalue), colour=sig)) +
  geom_point(size=1) +
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50')+
  guides(colour=guide_legend(title = NULL)) +
  scale_colour_manual(breaks = c("FALSE","TRUE"), values=c("grey33","red"), labels=c(expression("padj" >= "0.05"), "padj < 0.05"))+
  geom_text(aes(y=-log10(pvalue)+2, label=Gene), colour="black", size=2.5)



dev.off()



##Plot heatmap of vsd
pdf(file="Heatmap_of_VSD.pdf")
pheatmap(assay(vsd), show_rownames = FALSE, annotation_col = coldata, main = "Heatmap of transciptome expression (VSD)", show_colnames = FALSE)
dev.off()




## Sample to Sample distances
sampleDists <- dist(t(assay(vsd)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste( vsd$Tissue, vsd$Person, sep="-" )
colnames(sampleDistMatrix) <- NULL
cl <- colorRampPalette( rev(brewer.pal(9,"Blues"))) (255)
pdf(file="Sample_Distances_VSD.pdf")
pheatmap( sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col=cl, main = "Sample-Sample distances (VSD)")
dev.off()




## All data


## Subsetting Genes of intrest
## Extract list of genes which are sig DE in cortex (ie LFC > 0)
cortex <- na.omit(DiffExpression[ DiffExpression$`log2FoldChange (cortex over medulla)` > 0 & DiffExpression$padj < 0.05 ,])
## Extract list of genes which are sig DE in Meddulla (ie LFC < 0)
medula <- na.omit(DiffExpression[ DiffExpression$`log2FoldChange (cortex over medulla)` < 0 & DiffExpression$padj < 0.05 ,])

####### GSEA ############################################## 
## list of gene names and rank metiric (1/pvalue x abs(logfoldchange)logfoldchange)
GSEA <- na.omit(DiffExpression)
GSEA$gene <- row.names(GSEA)
GSEA$rank <- 1/GSEA$pvalue * abs( GSEA$log2FoldChange)/GSEA$log2FoldChange
GSEA <- GSEA[ , c("gene","rank")]
write.table(GSEA, file="GSEA.rnk", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)


## Pull out genes of intrest
genelist <- read.csv("genes_of_intrest.csv")
norm <- assay(vsd)
norm <- norm[,c(1,3,5,7,9,2,4,6,8,10)] ## reorder to give cortex then medulla
coldata2 <- coldata[c(1,3,5,7,9,2,4,6,8,10),]
gene1 <- norm[ row.names(norm) %in% genelist$V1,]
gene2 <- norm[ row.names(norm) %in% genelist$V2,]
gene3 <- norm[ row.names(norm) %in% genelist$V3,]
gene4 <- norm[ row.names(norm) %in% genelist$V4,]
gene5 <- norm[ row.names(norm) %in% genelist$V5,]
gene5 <- rbind( gene5, norm["FUT4",])
row.names(gene5)[10] <- "FUT4"
DE1 <- DiffExpression[ row.names(DiffExpression) %in% genelist$V1,]
DE1[is.na(DE1)] <- "FALSE"
DE2 <- DiffExpression[ row.names(DiffExpression) %in% genelist$V2,]
DE2[is.na(DE2)] <- "FALSE"
DE3 <- DiffExpression[ row.names(DiffExpression) %in% genelist$V3,]
DE3[is.na(DE3)] <- "FALSE"
DE4 <- DiffExpression[ row.names(DiffExpression) %in% genelist$V4,]
DE4[is.na(DE4)] <- "FALSE"
DE5 <- DiffExpression[ row.names(DiffExpression) %in% genelist$V5,]
DE5 <- rbind(DE5, DiffExpression["FUT4",])
DE5[is.na(DE5)] <- "FALSE"

## Make sumary for Gallager
Gal <- read.csv("genelist_gallagher.csv", header=FALSE)
gene6 <- norm[ row.names(norm) %in% Gal$V1,]
DE6 <- DiffExpression[ row.names(DiffExpression) %in% Gal$V1,]
gene6 <- merge( gene6,DE6, by =0 )
gene6 <- gene6[order(gene6$log2FoldChange, decreasing = TRUE),]
gene6 <- na.omit(gene6)
row.names(gene6) <- gene6$Row.names






## Make heatmap of genes of intrest
##Plot heatmap

## colour palate
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))
##Heatmap
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))
color.map2 <- function(sig) {if (sig=="TRUE") "red" else "black"}
sigcol <- unlist(lapply(DE1$sig, color.map2))

pdf(file="Heatmap_GS1.pdf")
heatmap.2(as.matrix(gene1), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          RowSideColors = sigcol
          )
legend("topright",
       legend=c("Cortex", "Medulla", expression("padj" >= "0.05"), "padj < 0.05"),
       col=c("grey","purple","black","red"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()



sigcol <- unlist(lapply(DE2$sig, color.map2))

pdf(file="Heatmap_GS2.pdf")
heatmap.2(as.matrix(gene2), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          RowSideColors = sigcol
)
legend("topright",
       legend=c("Cortex", "Medulla", expression("padj" >= "0.05"), "padj < 0.05"),
       col=c("grey","purple","black","red"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

sigcol <- unlist(lapply(DE3$sig, color.map2))

pdf(file="Heatmap_GS3.pdf")
heatmap.2(as.matrix(gene3), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          margin=c(10,10),
          ColSideColors=tissuecol,
          RowSideColors = sigcol
)
legend("topright",
       legend=c("Cortex", "Medulla", expression("padj" >= "0.05"), "padj < 0.05"),
       col=c("grey","purple","black","red"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()


sigcol <- unlist(lapply(DE4$sig, color.map2))

pdf(file="Heatmap_GS4.pdf")
heatmap.2(as.matrix(gene4), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          RowSideColors = sigcol
)
legend("topright",
       legend=c("Cortex", "Medulla", expression("padj" >= "0.05"), "padj < 0.05"),
       col=c("grey","purple","black","red"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

sigcol <- unlist(lapply(DE5$sig, color.map2))

pdf(file="Heatmap_GS5.pdf")
heatmap.2(as.matrix(gene5), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          RowSideColors = sigcol
)
legend("topright",
       legend=c("Cortex", "Medulla", expression("padj" >= "0.05"), "padj < 0.05"),
       col=c("grey","purple","black","red"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

sigcol <- unlist(lapply(gene6$sig, color.map2))
gene6 <- gene6[,c(2:11)]

pdf(file="Heatmap_Gal.pdf")
heatmap.2(as.matrix(gene6), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          RowSideColors = sigcol
)
legend("topright",
       legend=c("Cortex", "Medulla", expression("padj" >= "0.05"), "padj < 0.05"),
       col=c("grey","purple","black","red"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

## Data for Gal
# VSD
write.csv(gene6, file = "VSD_gal.csv")
# Normalized counts
normgal <- counts(dds, normalized = TRUE)
normgal <- normgal[row.names(normgal) %in% row.names(gene6),]
write.csv(normgal, file = "Norm_counts_gal.csv")
##Heat map of top and bottom 50 proteins ranked by DE


Resualts <- merge(as.data.frame(norm), DiffExpression, by="row.names")
row.names(Resualts) <- Resualts$Row.names
Resualts <- Resualts[,c(2:11,13)]

top <- Resualts[order(Resualts$log2FoldChange, decreasing = TRUE),]
top <- top[c(1:50),]
top <- top[,c(1:10)]
bottom <- Resualts[order(Resualts$log2FoldChange, decreasing = FALSE),]
bottom <- bottom[c(1:50),]
bottom <- bottom[,c(1:10)]

coldata2 <- coldata[c(1,3,5,7,9,2,4,6,8,10),]
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))

pdf(file="Cortex_DE_RNA.pdf")
heatmap.2(as.matrix(top), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

pdf(file="Medulla_DE_RNA.pdf")
heatmap.2(as.matrix(bottom), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()













coldata2 <- coldata[c(1,3,5,7,9,2,4,6,8,10),]
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))

pdf(file="Cortex_DE_Proteins.pdf")
heatmap.2(as.matrix(top), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

###########
##Comparision between RNASeq and Microarray
##Load in data from DE of microarray

Micro <- read.csv("~/Documents/LMB/Microarrary analysis/GSE3931_Kidney_Microarray/Diffrential_Expression_Kidney_Microarray.csv", row.names=1)
Micro2 <- read.csv("~/Documents/LMB/Microarrary analysis/GSE3931_Kidney_Microarray/Diffrential_Expression_Kidney_Microarray.csv", row.names=1)

## Make list of genes common to both microarray and seq

Micro <- Micro[row.names(Micro) %in% row.names(DiffExpression),]
Seq <- DiffExpressionNA[row.names(DiffExpressionNA) %in% row.names(Micro),]

comp <- merge( Micro, Seq, by="row.names")
comp <- comp[,c(1,2,6,9,13)]
row.names(comp) <- comp$Row.names
comp$Row.names <- NULL
colnames(comp) <- c("Micro_LFC", "Micro_padj", "Seq_LFC", "Seq_padj")
comp <- na.omit(comp)

sink("Micro_vs_rnaseq.txt")
cor.test(comp$Micro_LFC, comp$Seq_LFC)
sink()

comp$Msig <- apply(comp, 1, function(x)
{ if (as.numeric(x[2]) < 0.05){
  y<-10}
  else {
    y<-0}
  return(y)
})
comp$Ssig <- apply(comp, 1, function(x)
{ if (as.numeric(x[4]) < 0.05){
  y<-1}
  else{ 
    y<-0}
  return(y)
})
comp$Sig <- comp$Msig + comp$Ssig

pdf(file="Microarray_vs_RNASeq_Scatter.pdf")  
ggplot (data = comp, aes(x=Micro_LFC, y=Seq_LFC )) +
  geom_point(aes(color=as.factor(Sig))) +
  scale_colour_manual(breaks = c("0","1","10","11"), values=c("black", "green", "purple", "red"), labels=c(expression("Both padj" >= "0.05"), "RNASeq padj < 0.05", "Microarray padj < 0.05", "Both padj < 0.05"))+
  guides(color=guide_legend(title=NULL))+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Microarray (LFC)") +
  ylab("RNASeq (LFC)") +
  annotate("text", x=3, y=(-4), label="r= 0.7463968", hjust=0)+
  annotate("text", x=3, y=(-4.3), label="p< 2.2e-16", hjust=0)
dev.off()
##Venn diagram of overlap
pdf(file ="Venn_Microarray_vs_Seq.pdf")
draw.pairwise.venn(area1= nrow(Micro2), 
                   area2 = nrow(DiffExpressionNA), 
                   cross.area = nrow(comp), 
                   category = c( "Microarray", "RNASeq"),
                   lty=rep("blank",2),
                   fill = c("lightblue", "pink"),
                   cat.pos = c(0, 0), 
                   cat.dist = 0.01,
                   cat.fontface = 2)
dev.off()

## ven of significance for shared genes

A <- nrow(comp[comp$Msig > 0,]) ## number of significant DE in microarray
B <- nrow(comp[comp$Ssig > 0,]) ## number of significant DE in seq
c <- nrow(comp[comp$Sig == 11,]) ## Number significant for both 
pdf(file="Venn_Microarray_vs_seq_sig.pdf")
draw.pairwise.venn(area1= A, 
                  area2 = B, 
                  cross.area = c, 
                  category = c( "Microarray", "RNASeq"),
                  lty=rep("blank",2),
                  fill = c("lightblue", "pink"),
                  cat.pos = c(0, 0), 
                  cat.dist = 0.01,
                  cat.fontface = 2)
dev.off()


## Plot dispersion estimates

pdf(file="Dispersion_est.pdf")
plotDispEsts(dds)
dev.off()


#### genes which can be used to define if a sample is in the cortex or medulla

cvm <- na.omit(DiffExpression)
## filter for siginifance
cvm <- cvm[cvm$sig == "TRUE" ,]
## Combine with expression data
normcounts <- as.data.frame(counts(dds, normalized = TRUE))
cvm <- merge(cvm, normcounts, by=0)
cvm$cortex_Median <- rowMedians(as.matrix(cvm[,c(9,11,13,15,17)]))
cvm$medulla_Median <- rowMedians(as.matrix(cvm[,c(10,12,14,16,18)]))
## Filter for samples which have a median value >200 for one group and < 25 for the other
cvm <- cvm[cvm$cortex_Median > 200 & cvm$medulla_Median < 25 | cvm$medulla_Median > 200 & cvm$cortex_Median < 25,]
## Lable for tissue which the gene belogs to
cvm$tissue <- as.factor(apply(cvm,1, function(x) {if (as.numeric(x[19]) < 25)("Medulla") else ("Cortex")}))
## export gene list
tmp <- cvm[,c(19,20,22)]

## Filter for  |LFC| > 2.5
cvm <- cvm[ abs(cvm$log2FoldChange) > 2.5,] 
write.csv(tmp, file="~/Documents/Gene lists/Cortex_Vs_Medulla_genelist.csv")
## Plot filtered median counts against eachother
pdf(file="Cor_vs_Med_Log_Median_genes_filtered.pdf")
ggplot(data = cvm, aes(x=log10(cortex_Median), y=log10(medulla_Median)))+
  geom_point()+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Cortex - log10(Median Reads per gene)") +
  ylab("Medulla - log10(Median Reads per gene)")
dev.off()
 
pdf(file="Cor_vs_Med_Median_genes_filtered.pdf") 
ggplot(data = cvm, aes(x=cortex_Median, y=medulla_Median))+
  geom_point(size=7, shape=1)+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Cortex Median Reads per gene") +
  ylab("Medulla Median Reads per gene")
dev.off()


## For PCR validation could try top 5 in each group based on top median expression in each group
row.names(cvm) <- cvm$Row.names
Targetgenes <- c(row.names(cvm[order(cvm$cortex_Median, decreasing = TRUE),][c(1:5),]),row.names(cvm[order(cvm$medulla_Median, decreasing = TRUE),][c(1:5),]))
Targetgenes <- DiffExpression[row.names(DiffExpression) %in% Targetgenes,]
write.csv(Targetgenes, file="CvM_qPCR_TARGETS.csv")

## re plot graph with dots coloured for presence in target gene set

cvm$target <- as.factor(row.names(cvm) %in% row.names(Targetgenes))

pdf(file="Cor_vs_Med_Median_genes_filtered_coloured.pdf") 
ggplot(data = cvm, aes(x=cortex_Median, y=medulla_Median, color=target))+
  geom_point()+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Cortex - Median Reads per gene") +
  ylab("Medulla - Median Reads per gene")+
  scale_colour_manual(breaks = c("FALSE","TRUE"), values=c("grey33","red"))+
  theme(legend.position="none")
dev.off()


## Heatmap of filtered genes

gene7 <- norm[ row.names(norm) %in% cvm$Row.names,]
DE7 <- DiffExpression[ row.names(DiffExpression) %in% cvm$Row.names,]
gene7 <- merge( gene7,DE7, by =0 )
gene7 <- gene7[order(gene7$log2FoldChange, decreasing = TRUE),]
gene7 <- na.omit(gene7)
row.names(gene7) <- gene7$Row.names

gene7 <- gene7[,c(2:11)]
mybreaks = mybreaks <- c(seq(5,7,length=50),seq(7.1,10,length=50))

pdf(file="Heatmap_cvm.pdf")
heatmap.2(as.matrix(gene7), 
          trace="none", 
          tracecol="black", 
          col=hmcol,
          breaks = mybreaks,
          key.title ="", 
          key.ylab="", 
          key.xlab="VSD expression",
          key=TRUE, 
          scale="none", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

##### GSEA graphs
## Hallmarks
Hall <- rbind( read.csv("~/Documents/LMB/BTRU/Kidney/Non perfused study/RNASeq_April_2016/GSEA/Hallmarks/gsea_report_neg.csv"),read.csv("~/Documents/LMB/BTRU/Kidney/Non perfused study/RNASeq_April_2016/GSEA/Hallmarks/gsea_report_pos.csv"))
Hall <- Hall[Hall$FDR.q.val < 0.05,]
tmp <- substr(Hall$NAME, 10,100)
tmp <- gsub("_", " ", tmp)
row.names(Hall) <- tmp
Hall <- Hall[,c(6,8)]
Hall$Name <- factor(row.names(Hall),levels = unique(row.names(Hall)) )
Hall <- Hall[order(Hall$NES, decreasing = FALSE),]

## Plot bar graphs
pdf(file="GSEA_Hallmarks.pdf", height = 4)
ggplot( data=Hall, aes(x=Name, y=NES, fill=FDR.q.val)) +
  geom_bar(stat="identity") +
  coord_flip()+
  geom_hline(yintercept=0, colour="grey50") +
  scale_fill_gradient(low="orangered2", high="black") +
  labs(y="NES", x="Hallmarks")
dev.off()


Hallmarks <- rbind( read.csv("~/Documents/LMB/BTRU/Kidney/Non perfused study/RNASeq_April_2016/GSEA/Hallmarks/gsea_report_neg.csv"),read.csv("~/Documents/LMB/BTRU/Kidney/Non perfused study/RNASeq_April_2016/GSEA/Hallmarks/gsea_report_pos.csv"))
Kegg <- rbind( read.csv("~/Documents/LMB/BTRU/Kidney/Non perfused study/RNASeq_April_2016/GSEA/Kegg/gsea_report_neg.csv"),read.csv("~/Documents/LMB/BTRU/Kidney/Non perfused study/RNASeq_April_2016/GSEA/Kegg/gsea_report_pos.csv"))
tmp <- substr(Hallmarks$NAME, 10,100)
tmp <- gsub("_", " ", tmp)
row.names(Hallmarks) <- tmp
Hallmarks <- Hallmarks[,c(6,8)]
tmp <- substr( Kegg$NAME, 6,100)
tmp <- gsub("_", " ", tmp)
row.names(Kegg) <- tmp
Kegg <- Kegg[,c(6,8)]
## Remove any sets which are non signifant (FDR > 0.05)
Hallmarks <- Hallmarks [Hallmarks$FDR.q.val< 0.05,]
Hallmarks <- Hallmarks[order(Hallmarks$NES, decreasing = TRUE),]
Hallmarks$Hallmark <- factor(row.names(Hallmarks), levels = unique(row.names(Hallmarks)))
Kegg <- Kegg [Kegg$FDR.q.val< 0.05,]
Kegg <- Kegg[order(Kegg$NES, decreasing = TRUE),]
Kegg$Kegg <- factor(row.names(Kegg), levels = unique(row.names(Kegg)))
## Plot bar graphs
pdf(file="GSEA_Hallmarks.pdf", useDingbats = FALSE)
ggplot( data=Hallmarks, aes(x=Hallmark, y=NES, fill=FDR.q.val)) +
  geom_bar(stat="identity") +
  coord_flip()+
  geom_hline(yintercept=0, colour="grey50") +
  scale_fill_gradient(low="orangered2", high="black") +
  labs(y="NES", x="Hallmarks")
dev.off()
pdf(file="GSEA_Kegg.pdf",width = 12, useDingbats = FALSE)
ggplot( data=Kegg, aes(x=Kegg, y=NES, fill=FDR.q.val)) +
  geom_bar(stat="identity") +
  coord_flip()+
  geom_hline(yintercept=0, colour="grey50") +
  scale_fill_gradient(low="orangered2", high="black") +
  labs(y="NES", x="Kegg")
dev.off()


## plot counts for menna

d <- plotCounts(dds, gene = "NALCN", intgroup = "Tissue", returnData = TRUE)
d <- merge (d, coldata, by=0)
pdf("NALCN gene expression.pdf")
ggplot(data = d, aes(x=Tissue.x , y=count, group = Person))+
  geom_point(size = 4)+
  geom_line( colour="black", linetype = 3, size=1) +
  xlab("Tissue") +
  ylab("Count")+
  labs(title = "NALCN - RNASeq data")+
  geom_hline(yintercept=0, color='grey50')
dev.off

d <- d[! d$Person == "A",] ## remove the sample with contamination
pdf("NALCN gene expression.pdf")
ggplot(data = d, aes(x=Tissue.x , y=count, group = Person))+
  geom_point(size = 4)+
  geom_line( colour="black", linetype = 3, size=1) +
  xlab("Tissue") +
  ylab("Count")+
  labs(title = "NALCN - RNASeq data")+
  geom_hline(yintercept=0, color='grey50')
dev.off()


## Heatmap of Heamoglobin genes
Hb <- norm[grep("HB", substr( row.names(norm), 1, 2)),] ## Pull out any gene where the first two letters of the gene name are HB
tmp <- c("HBEGF","HBS1L","HBP1")
Hb <- Hb[! row.names(Hb) %in% tmp,]

## colour palate
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))
##Heatmap
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))
## Remove contaminated samples and reorder rows
Hb <- Hb[,c(-1,-6)]
Hb <- Hb[c(7,6,1,3,4,2,5),]
tissuecol <- tissuecol[c(-1,-6)]
pdf(file="Heatmap_Hb.pdf")
heatmap.2(as.matrix(Hb), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()


## Heatmap of genes for menna 

tmp <- c("TGFB1","AREG","ADAM17")
Hb <- norm[row.names(norm) %in% tmp,]

## colour palate
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))
##Heatmap
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))
## Remove contaminated samples and reorder rows
Hb <- Hb[,c(-1,-6)]
tissuecol <- tissuecol[c(-1,-6)]
pdf(file="Heatmap_AREG_Kidney_RNASEQ.pdf")
heatmap.2(as.matrix(Hb), 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          cexRow = 3,
          margins = c(6,12)
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7
)
dev.off()

## Boxplot of HBGEF
tmp <- as.data.frame(norm)
tmp <- tmp["HBEGF",]
tmp <- as.data.frame(t(tmp))
tmp$tissue <- c(rep("Cortex",5),rep("Medulla", 5))

pdf(file="Boxplot_HBEGF_RNASeq.pdf", useDingbats = FALSE)
boxplot(tmp$HBEGF ~ tmp$tissue)
text(x=2, y=7.4, labels = "padj = 0.8060499")
dev.off()
save.image()


## Make heatmap of top 10 genes for each directions
tmp <- na.omit(DiffExpression[DiffExpression$padj < 0.05,])
tmp <- tmp[order(tmp$log2FoldChange),]
tmp <- tmp[c(1:10,(nrow(tmp)-9):nrow(tmp)),]
tmp <- merge (tmp, norm, by = 0)
row.names(tmp) <- tmp$Row.names
tmp <- tmp[order(tmp$log2FoldChange),]
tmp <- as.matrix(tmp[,c(9:18)])
## colour palate
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))
##Heatmap
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))


pdf(file="Heatmap_top10_eachway.pdf")
heatmap.2(tmp, 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          margins = c(3,7)
)
legend("topright",
       legend=c("Cortex", "Medulla"),
       col=c("grey","purple"),
       lty=1,
       lwd=5,
       cex=.7,
       bg = "white"
)
dev.off()


## NK Chemokines

NK <- read.csv("NK_CHEMOKINES.CSV", header = FALSE)
tmp <- na.omit(DiffExpression[row.names(DiffExpression) %in% NK$V1,])
tmp <- merge (tmp, norm, by = 0)
row.names(tmp) <- tmp$Row.names
tmp <- tmp[order(tmp$log2FoldChange),]
sigcol <- unlist(lapply(tmp$sig, color.map2))
tmp <- as.matrix(tmp[,c(9:18)])
## colour palate
hmcol <- c(colorRampPalette(c("blue","white","red"))(99))
##Heatmap
color.map <- function(Tissue) {if (Tissue=="cortex") "grey" else "purple"}
tissuecol <- unlist(lapply(coldata2$Tissue, color.map))


pdf(file="Heatmap_NK.pdf")
heatmap.2(tmp, 
          trace="none", 
          tracecol="black", 
          col=hmcol,  
          key.title ="", 
          key.ylab="", 
          key.xlab="Row Z score",
          key=TRUE, 
          scale="row", 
          Colv=FALSE, Rowv=FALSE, dendrogram="none", 
          labCol = FALSE,
          ColSideColors=tissuecol,
          RowSideColors = sigcol,
          margins = c(3,7)
)
legend("topright",
       legend=c("Cortex", "Medulla", "Sig", "not Sig"),
       col=c("grey","purple", "red", "black"),
       lty=1,
       lwd=5,
       cex=.7,
       bg = "white"
)
dev.off()

## Analysis of key genes for the validation of the single cell paper Feb 2018

GOI <- c("CXCL12", "CXCL14")

tmp <- norm[row.names(norm) %in% GOI,]
colnames(tmp) <- c(rep("cortex", 5), rep("medulla",5))

pdf("~/Dropbox/Figures for scRNAseq kidney paper/Heatmap_CXCL12_CXCL14_Kidney_ourdata.pdf")
pheatmap(tmp, cluster_cols = FALSE, scale = "row" )
dev.off()

GOI <- c("CXCL5", "CXCL2", "CXCL8", "CXCL3", "CXCL7")
tmp <- norm[row.names(norm) %in% GOI,]
colnames(tmp) <- c(rep("cortex", 5), rep("medulla",5))

pdf("~/Dropbox/Figures for scRNAseq kidney paper/Heatmap_CXCR2_ligands_Kidney_ourdata.pdf")
pheatmap(tmp, cluster_cols = FALSE, scale = "row" )
dev.off()



GOI <- c("REG3G", "REG3A", "SAA1", "DEFA5", "DEFB1", "LCN2", "UMOD", "CAMP", "PTX3")
tmp <- norm[row.names(norm) %in% GOI,]
colnames(tmp) <- c(rep("cortex", 5), rep("medulla",5))

pdf("~/Dropbox/Figures for scRNAseq kidney paper/Heatmap_AMP_Kidney_ourdata.pdf")
pheatmap(tmp, cluster_cols = FALSE, scale = "row")
dev.off()

## All cytokines
cytokines <- read.csv("~/Documents/Gene lists/Human_Cytokines_V3_no_Chemokines.csv")
tmp <- norm[row.names(norm) %in% cytokines$Cytokines,]
colnames(tmp) <- c(rep("cortex", 5), rep("medulla",5))

pdf("~/Dropbox/Figures for scRNAseq kidney paper/Heatmaps from John/Heatmap_Cytokines_ourdata.pdf", onefile = FALSE)
pheatmap(tmp, cluster_cols = FALSE, scale = "row")
dev.off()

## all chemokines and chemokine receptors
tmp <- readRDS("~/Documents/Gene lists/chemokines etc genesets.RDS")
chemokines <- c(tmp$Chemokines, tmp$ChemokineR)
tmp <- norm[row.names(norm) %in% chemokines,]
colnames(tmp) <- c(rep("cortex", 5), rep("medulla",5))

pdf("~/Dropbox/Figures for scRNAseq kidney paper/Heatmaps from John/Heatmap_chemokines_ourdata.pdf", onefile = FALSE)
pheatmap(tmp, cluster_cols = FALSE, scale = "row")
dev.off()



############### Use single cell data to look at expression of HVG for Podocytes and colecting duct cells

HVG <- readRDS("nephron cluster markers lists.RDS")
tmp <- levels(HVG$cluster)
tmp2 <- coldata[,-3]
for (i in tmp){
  a <- as.character(HVG$gene.symbol[HVG$cluster == i])
  a <- norm[row.names(norm) %in% a,]
  pdf(file = paste("Comparision_to_HVG/Heatmap_",i,".pdf", sep=""), onefile = FALSE, useDingbats = FALSE)
  print(pheatmap(a, scale="row", show_colnames = FALSE, annotation_col = tmp2))
  dev.off()
}


tmp <- as.character(HVG$gene.symbol[HVG$cluster == "Podocytes"])

tmp <- norm[row.names(norm) %in% tmp,]
tmp2 <- coldata[,-3]

pheatmap(tmp, scale="row", show_colnames = FALSE, annotation_col = tmp2)

######## CvM - make data frame of VSD which can be used for making a model

tmp <- as.data.frame(assay(vsd))
tmp <- tmp[row.names(tmp) %in% row.names(cvm),]
tmp <- as.data.frame(t(tmp))
tmp$Tissue <- coldata$Tissue
## remove the sample which we know is not pure (persone A medulla sample)
tmp <- tmp[-2,]

write.csv(tmp, file = "~/Documents/Gene lists/CvM_training_set.csv")
