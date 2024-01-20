#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

setwd(".")
library("DESeq2") # V.1.34.0
library(ggdendro) # V0.1.23

set.seed(0)

# plot dendrogram for paired kidney samples

# load count data
load("gene_counts_table.RData")
countTable <- counts_tab$counts
tmp <- colnames(countTable)
tmp <- strsplit(tmp, "_")
tmp <- lapply(tmp, function(x) {x[1]})
colnames(countTable) <- tmp
countTable <- countTable[,order(as.numeric(colnames(countTable)))]

## Sample information
Person <- as.factor(c("5", "3", "1", "4", "2", "3", "2", "3", "5", "1", "5", "2", "4", "4", "3","5","2", "1", "4", "1"))
Condition <- factor(c("pre", "pre", "pre", "post", "pre", "post", "post", "post", "pre","post", "post", "pre","pre","post", "pre","post", "post", "post", "pre","pre"), levels = c("pre", "post"))
Filter <- factor(c("Y", "Y","N","N","N","N","N","Y","N","Y","N","Y","Y","Y","N", "Y","Y","N","N", "Y"), levels = c("N","Y"))
coldata <- data.frame(Person, Condition, Filter, row.names = colnames(countTable))
coldata$pairs <- as.factor(paste(coldata$Person, coldata$Filter, sep="."))
coldata$Sex <- as.factor(apply(coldata,1, function(x) {if (x[1] == "5" | x[1] == "3"| x[1] == "4")("M") else ("F")}))

# DE analysis
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = ~Person + Filter + Condition)
dds <- dds [ rowSums(counts(dds)) >1 ,]
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
resualtsPerf <- as.data.frame(res)
resualtsPerf$sig <- as.factor(resualtsPerf$padj < 0.05)

# Variance stabilised transformation
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
tmp <- assay(vsd)
tmp <- tmp[,coldata$Condition == "pre"]
coldata <- coldata[coldata$Condition == "pre",] 
lab <- c("5A",
         "3A",
         "1B",
         "2B",
         "5B",
         "2A",
         "4A",
         "3B",
         "4B",
         "1A")
colnames(tmp) <- lab

# create dendrogram
sampleTree <- hclust(dist(t(tmp)))
dend <- as.dendrogram(sampleTree)
dend_data <- dendro_data(dend, type = "rectangle")
head(dend_data$segments)
head(dend_data$labels)
ggdendrogram(sampleTree)