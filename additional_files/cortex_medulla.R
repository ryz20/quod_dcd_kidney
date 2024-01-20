#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

library(tidyverse) # V1.3.2
library(DESeq2) # V1.34.0

set.seed(1)

# Initial analysis by J R Ferdinand commencing May 2016 using DESeq2 for DE.
# Results and code can be found in 'Additional Files'

# Import files
RNAseq_DE <- read.csv(file = "CvM_RNAseq_DE.csv") # RNAseq DE results
row.names(RNAseq_DE) <- RNAseq_DE$X
RNAseq_DE$X <- NULL
colnames(RNAseq_DE)[2] <- "log2FoldChange"
load("gene_counts_table.RData") # Raw counts
countTable <- counts_tab$counts
colnames(countTable) <- c(1,10,2,3,4,5,6,7,8,9)
countTable <- countTable [, c(1,3,4,5,6,7,8,9,10,2)]
coldata <- read.csv("RNAseqdesign.csv", row.names = 1) # Metadata
coldata$Sample <- 1:10

# Create cortex/medulla specific genelist
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = ~Person + Tissue)
dds <- dds [ rowSums(counts(dds)) >1 ,]
dds$Tissue <- relevel(dds$Tissue, ref="medulla")
dds <- DESeq(dds, betaPrior = T)
normcounts <- as.data.frame(counts(dds, normalized = TRUE))
# Select CvM genelist
RNAseq_DE$sig <- as.factor(RNAseq_DE$padj < 0.05)
cvm <- na.omit(RNAseq_DE)
cvm <- cvm[cvm$sig == "TRUE" ,]
normcounts <- normcounts[match(row.names(cvm), row.names(normcounts)),]
cvm <- merge(cvm, normcounts, by=0)
# Calculate cortex/medulla medians (due to outlier)
cvm$cortex_Median <- rowMedians(as.matrix(cvm[,c(9,11,13,15,17)]))
cvm$medulla_Median <- rowMedians(as.matrix(cvm[,c(10,12,14,16,18)]))
# Filter for samples which have a median value >200 for one group and < 25 for the other
cvm <- cvm[cvm$cortex_Median > 200 & cvm$medulla_Median < 25 | cvm$medulla_Median > 200 & cvm$cortex_Median < 25,]
# Lable for tissue which the gene belogs to
cvm$tissue <- as.factor(apply(cvm,1, function(x) {if (as.numeric(x[19]) < 25)("Medulla") else ("Cortex")}))
# Filter for  |LFC| > 2.5
cvm <- cvm[abs(cvm$log2FoldChange) > 2.5,] 
row.names(cvm) <- cvm$Row.names
cvm$Row.names <- NULL
# export gene list
write.csv(cvm[,c(18,19,20)], file="Cortex_Vs_Medulla_genelist.csv")


# Data visualization
table(RNAseq_DE$sig) # 5195 total significant genes
table(RNAseq_DE[RNAseq_DE$log2FoldChange > 0,]$sig) # 2913 significant up genes
table(RNAseq_DE[RNAseq_DE$log2FoldChange < 0,]$sig) # 2282 significant up genes
# Volcano plot
res.volcano <- na.omit(as.data.frame(RNAseq_DE))
sigGenes <- res.volcano %>% 
  filter(padj<0.05) %>% 
  rownames_to_column(var = "gene") %>% 
  pull("gene")
write.csv(sigGenes, file = "RNAseq_DE.csv")
res.volcano <- res.volcano %>% 
  mutate(sig = ifelse(padj < 0.05, "padj<0.05", "padj>0.05"))
lab <- res.volcano$padj < 0.05
res.volcano$lab <- row.names(res.volcano)
res.volcano$lab[!lab] <- NA
pdf(file="RNAseq_volcano.pdf", useDingbats = FALSE)
ggplot(data = res.volcano, aes(x=log2FoldChange, y=-log10(pvalue), colour=sig, label=lab)) +
  geom_point(size=2) +
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50')+
  guides(colour=guide_legend(title = NULL)) +
  labs(x= "Log2 Fold Change", y="-log10 (p value)")
dev.off() 
# Heatmap of top DE genes
res.heatmap <- as.data.frame(RNAseq_DE)
# Filter DE genes, select top n (e.g. 100) and pull names
res.heatmap <- rownames_to_column(res.heatmap, var = "gene")
sigGenes <- res.heatmap %>% 
  filter(padj < 0.05) %>% 
  arrange(desc(log2FoldChange)) 
sigGenes <- sigGenes[c(1:10, (nrow(sigGenes)-9):nrow(sigGenes)),]
sigGenes <- sigGenes %>% 
  pull("gene")
## Extract transformed counts for these genes
plotDat <- vst(dds, blind = FALSE)[sigGenes,] %>% 
  assay()
# Set annotations
anno_info <- as.data.frame(colData(dds)[,c("Sample", "Tissue")])
# Change order of samples
anno_info <- anno_info[order(anno_info$Tissue),]
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
threshold <- 2
plotDat2 <-  apply(plotDat2, c(1,2), trim)
# plot heatmap
pdf(file="RNAseq_heatmap.pdf", useDingbats = FALSE)
pheatmap::pheatmap(mat = plotDat, 
         scale="row", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = anno_info)
pheatmap::pheatmap(mat = plotDat2, 
         scale="none", 
         cluster_rows = FALSE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE,
         show_colnames = TRUE,
         annotation_col = anno_info)
dev.off()


# Correlation with Mass Spec LFC
# load mass spec DE
MassSpec_DE <- read.csv(file = "CvM_MassSpec_DE.csv") # 4288 proteins
RNA_Counts <- read.csv(file = "VSD_counts.csv", row.names = 1) # variance stabilised transformed values
RNA <- merge(RNAseq_DE, RNA_Counts, by="row.names")
rm(RNA_Counts)
row.names(RNA) <- RNA$Row.names
RNA$Row.names <- NULL
# Filter for genes only found in both
MSfilt <- MassSpec_DE[MassSpec_DE$Gene %in% row.names(RNA),]
RNAfilt <- RNA[ row.names(RNA) %in% MassSpec_DE$Gene ,]
# 4 of the genes are duplicated in MSfilt, removed these
MSfilt <- MSfilt[! duplicated(MSfilt$Gene), ]
# Make gene names row names
row.names(MSfilt) <- MSfilt$Gene
# Correlation of DE between samples
comb <- merge(MSfilt, RNAfilt, by="row.names" )
comb <- comb[, c(2:4,6,18,22)]
colnames(comb) <- c("Gene", "Protein", "Prot.LFC","Prot.padj", "RNA.LFC", "RNA.padj")
comb <- na.omit(comb)
# Add Columns to enable colouring by significance
comb$Msig <- apply(comb, 1, function(x)
{ if (as.numeric(x[4]) < 0.05){
  y<- 1}
  else {
    y <- 0}
  return(y)
})
comb$Rsig <- apply(comb, 1, function(x)
{ if (as.numeric(x[6]) < 0.05){
  y<- 10}
  else {
    y <- 0}
  return(y)
})
comb$sig <- comb$Msig + comb$Rsig
# calculate pearson's correlation
sink(file = "RNASeq_vs_Mass_Spec_LFC.txt")
cor.test(comb$Prot.LFC, comb$RNA.LFC, method = "pearson")
sink()
# plot outliers
model <- lm(Prot.LFC ~ RNA.LFC, data = comb)
summary(model)
rstandard(model)
comb$residuals <- rstandard(model)
comb$outlier <- abs(comb$residual) >2
comb$outlier.name <- ifelse(comb$outlier == TRUE, comb$Gene, NA)
outliers <- na.omit(comb$outlier.name)
write.table(outliers, file = "RNASeq_vs_Mass_Spec_LFC_outliers.tsv", quote = FALSE, sep = '\t', row.names = F, col.names = F)
# plot correlation between RNA and MS LFC
pdf(file="RNASeq_vs_Mass_Spec_LFC.pdf")
ggplot( data=comb, aes(x=Prot.LFC, y=RNA.LFC)) +
  geom_point(aes(color=as.factor(sig))) +
  scale_colour_manual(breaks = c("0","1","10","11"), values=c("black", "green", "purple", "red"), labels=c(expression("Both padj" >= "0.05"), "Protein padj < 0.05", "RNA padj < 0.05", "Both padj < 0.05"))+
  guides(color=guide_legend(title=NULL))+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Protein (LFC)") +
  ylab("RNA (LFC)") +
  annotate("text", x=2, y=(-3), label="r= 0.6832", hjust=0)+
  annotate("text", x=2, y=(-3.2), label="p= <2.2e-16", hjust=0)
ggplot(data=comb, aes(x=Prot.LFC, y=RNA.LFC, label = outlier.name)) +
  geom_point(aes(color=as.factor(sig))) +
  scale_colour_manual(breaks = c("0","1","10","11"), values=c("black", "green", "purple", "red"), labels=c(expression("Both padj" >= "0.05"), "Protein padj < 0.05", "RNA padj < 0.05", "Both padj < 0.05"))+
  guides(color=guide_legend(title=NULL))+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Protein (LFC)") +
  ylab("RNA (LFC)") +
  annotate("text", x=2, y=(-3), label="r= 0.6832", hjust=0)+
  annotate("text", x=2, y=(-3.2), label="p= <2.2e-16", hjust=0)+
  ggrepel::geom_text_repel(size = 2.5, box.padding = 0.1, lebel.padding = 0.1)
dev.off()

### GSEA Comparison with Mass Spec
neg <- read.csv("CvM_MassSpec_Hallmarks/gsea_report_neg.csv")
pos <- read.csv("CvM_MassSpec_Hallmarks/gsea_report_pos.csv")
HallMS<- rbind(neg, pos)
neg <- read.csv("CvM_RNASeq_Hallmarks/gsea_report_neg.csv")
pos <- read.csv("CvM_RNASeq_Hallmarks/gsea_report_pos.csv")
HallRNA<- rbind(neg, pos)
rm(pos)
rm(neg)
Hall<- merge(HallRNA, HallMS, by="NAME")
Hall <- Hall[,c(1,6,8,16,18)]
colnames(Hall) <- c("NAME", "NES_RNA", "FDR.q.val_RNA", "NES_Prot", "FDR.q.val_Prot")
hallmarks <- substr(Hall$NAME, 10, 60)
hallmarks <- gsub("_", " ", hallmarks)
row.names(Hall) <- hallmarks
Hall <- Hall[,(-1)]
Hall$Rsig <- apply(Hall, 1, function(x)
{ if (as.numeric(x[2]) < 0.05){
  y<-10}
  else {
    y<-0}
  return(y)
})
Hall$Msig <- apply(Hall, 1, function(x)
{ if (as.numeric(x[4]) < 0.05){
  y<-1}
  else{ 
    y<-0}
  return(y)
})
Hall$Sig <- Hall$Rsig + Hall$Msig
# calculate and plot hallmark NES correlations
sink(file="RNASeq_vs_Mass_Spec_hallmarks_correlation.txt")
cor.test(Hall$NES_Prot,  Hall$NES_RNA)
sink()
pdf(file="RNASeq_vs_Mass_Spec_hallmarks_scatter.pdf")  
ggplot (data = Hall, aes(x=NES_Prot, y=NES_RNA )) +
  geom_point(aes(color=as.factor(Sig))) +
  scale_colour_manual(breaks = c("0","1","10","11"), values=c("black", "green", "purple", "red"), labels=c(expression("Both padj" >= "0.05"), "Protein padj < 0.05", "RNA padj < 0.05", "Both padj < 0.05"))+
  guides(color=guide_legend(title=NULL))+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Protein NES") +
  ylab("RNA NES") +
  annotate("text", x=3, y=(-2.5), label="r= 0.5609622", hjust=0)+
  annotate("text", x=3, y=(-2.7), label="p= 5.01e-05", hjust=0)
dev.off()
# plot significantly enriched pathways in both RNA + MS
tmp <- Hall[Hall$Sig == 11 & (Hall$NES_RNA * Hall$NES_Prot) > 0,]
tmp1 <- tmp[,c(1,2)]
tmp1$omic <- "RNA" 
tmp1$Name <- row.names(tmp1)
colnames(tmp1) <- c("NES", "FDR", "Omic", "Name")
tmp2 <- tmp[,c(3,4)]
tmp2$omic <- "Protein" 
tmp2$Name <- row.names(tmp2)
colnames(tmp2) <- c("NES", "FDR", "Omic", "Name")
tmp<- rbind(tmp1, tmp2)
rm(tmp1)
rm(tmp2)   
pdf(file="RNASeq_vs_Mass_Spec_hallmarks_sig.pdf")
ggplot( data=tmp, aes(x=reorder(Name, NES), y=NES, fill=Omic), group=Omic) +
  geom_bar(stat="identity", position="dodge") +
  coord_flip()+
  geom_hline(yintercept=0, colour="grey50")+
  xlab("Hallmarks")+
  annotate("text", x="E2F TARGETS", y=(3), label="all FDR q value < 0.05", vjust=2)
dev.off()


# plot GSEA for RNA-seq alone
MasterHallmarks <- data.frame(matrix(vector(),0,4,dimnames = list(c(), c("NES", "FDR.q.val", "Comparision", "pathway"))),stringsAsFactors = FALSE)
tmp <- list.files(paste("CvM_RNAseq_Hallmarks/", sep=""))
tmp <- grep("gsea_report", tmp, value = TRUE)
tmp <- grep("csv", tmp, value = TRUE)
tmp <- paste("CvM_RNAseq_Hallmarks/",tmp, sep="")
tmp1 <- read.table(tmp[1], fill = TRUE, sep=",", header=TRUE)
tmp2 <- read.table(tmp[2], fill = TRUE, sep=",", header=TRUE)
Hallmarks <- rbind(tmp1, tmp2)
tmp <- substr(Hallmarks$NAME, 10,100)
tmp <- gsub("_", " ", tmp)
row.names(Hallmarks) <- tmp
Hallmarks <- Hallmarks[,c(5,6,8)]
Hallmarks$pathway <- row.names(Hallmarks)
MasterHallmarks <- rbind(MasterHallmarks, Hallmarks)
## Remove any sets which are non signifant (FDR > 0.05)
Hallmarks <- Hallmarks [Hallmarks$FDR.q.val< 0.05,]
Hallmarks <- Hallmarks[order(Hallmarks$NES, decreasing = FALSE),]
Hallmarks$Hallmark <- factor(row.names(Hallmarks), levels = unique(row.names(Hallmarks)))
Hallmarks <- droplevels(Hallmarks)
Hallmarks$dir <- factor(Hallmarks$NES > 0, levels = c("FALSE", "TRUE"))
if (all(Hallmarks$NES > 0)) {
  ## If all Hallmarks are positve
  pdf(file=paste("Hallmarks_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
  print(ggplot( data=Hallmarks, aes(x=Hallmark, y=NES,  size=-log10(FDR.q.val+0.00001))) +
          geom_point(colour="red") +
          geom_hline(yintercept=0, colour="grey50") +
          labs(y="NES", x="Hallmarks", size="-log10(FDR+0.0001)", colour="")+
          theme_bw()+
          coord_flip())
  dev.off()
}  else if (all(Hallmarks$NES < 0)){
  ## If all hallmarks are negative
  pdf(file=paste("Hallmarks_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
  print(ggplot( data=Hallmarks, aes(x=Hallmark, y=NES, size=-log10(FDR.q.val+0.00001))) +
          geom_point(colour="dodgerblue2") +
          geom_hline(yintercept=0, colour="grey50") +
          labs(y="NES", x="Hallmarks", size="-log10(FDR+0.0001)", colour="")+
          theme_bw()+
          coord_flip())
  dev.off()
}  else {
  ## Plot bar graphs
  pdf(file=paste("Hallmarks_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
  print(ggplot( data=Hallmarks, aes(x=Hallmark, y=NES, colour=dir, size=-log10(FDR.q.val+0.00001))) +
          geom_point() +
          geom_hline(yintercept=0, colour="grey50") +
          labs(y="NES", x="Hallmarks", size="-log10(FDR+0.0001)", colour="")+
          scale_colour_manual(breaks = c("FALSE","TRUE"), values=c("dodgerblue2","red"), labels=c("Downregulated\nPathway", "Upregulated\nPathway"))+
          theme_bw()+
          coord_flip())
  dev.off()
}


# Comparision between RNASeq and Microarray
# Load in data from DE of microarray (Higgins et al. 2004)
Micro_DE <- read.csv("Higgins_microarray_DE.csv", row.names=1)
library(ggVennDiagram) # V1.2.2
# list of RNAseq DEGs
tmp <- RNAseq_DE %>% 
  filter(sig == "TRUE")
sigRNAseq <- row.names(tmp)
# list of micro DEGs
tmp <- Micro_DE %>% 
  filter(adj.P.Val < 0.05)
tmp <- na.omit(tmp)
sigMicro <- row.names(tmp)
x <- list(RNA = sigRNAseq, Microarray = sigMicro)
p1 <- ggVennDiagram(x, label = "count", label_alpha = 0) +
  scale_fill_gradient (low = "#F4FAFE", high = "4981BF") +
  scale_color_manual(values = c("RNA" = "yellow","B" ="steelblue")) +
  theme(legend.position="none") +
  theme(text = element_text(size = 8,family = "sans")) 
# correlation
Micro_DE2 <- Micro_DE[row.names(Micro_DE) %in% row.names(RNAseq_DE),]
RNAseq_DE2 <- RNAseq_DE[row.names(RNAseq_DE) %in% row.names(Micro_DE),]
comp <- merge(Micro_DE2, RNAseq_DE2, by="row.names")
comp <- comp[,c(1,2,6,9,13)]
row.names(comp) <- comp$Row.names
comp$Row.names <- NULL
colnames(comp) <- c("Micro_LFC", "Micro_padj", "Seq_LFC", "Seq_padj")
comp <- na.omit(comp)
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
comp$BothSig <- comp$Msig + comp$Ssig
comp$sig <- ifelse(comp$BothSig == 11, "both sig",
                   ifelse(comp$BothSig == 0, "neither sig", "one sig"))
# statistics for annotation
sink("stats.txt")
cor.test(comp$Micro_LFC, comp$Seq_LFC)
sink()
# plot scatter
p1 <- ggplot (data = comp, aes(x=Micro_LFC, y=Seq_LFC )) +
  geom_point(aes(color=as.factor(sig)), shape = 1) +
  scale_colour_manual(values=c("#FF1F5B", "black", "#009ADE"), labels=c("Both", "None", "One"))+
  guides(color=guide_legend(title=NULL))+
  geom_hline(yintercept=0, color='grey50') +
  geom_vline(xintercept = 0, color='grey50') +
  xlab("Microarray LFC") +
  ylab("RNASeq LFC") +
  theme_bw() +
  theme(legend.position="none") +
  theme(text = element_text(size = 8,family = "sans")) 
