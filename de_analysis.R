#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

library(tidyverse) # V1.3.2
library(DESeq2) # V1.34.0

set.seed(1255342)

# DE analysis for n=271 samples

# Import metadata and count data
coldata <- read.csv("Coldata.csv")
row.names(coldata) <- coldata$X
coldata$X <- NULL
countTable <- read.csv("CountTable.csv")
row.names(countTable) <- countTable$X
countTable$X <- NULL
coldata[,'Sample_Set']<-factor(coldata[,'Sample_Set'])
coldata[,'Cluster']<-factor(coldata[,'Cluster'])
coldata <- droplevels(coldata)


# Differential expression analysis for each outcome accounting for missing data
# 1. DGF
coldata_DGF <- coldata %>% 
  mutate(DGF_bin = cut(DGF_days, breaks = c(-0.5, 1.5, Inf), 
                          labels = c("noDGF", "DGF"))) %>% 
  filter(!is.na(DGF_bin))
countTable_DGF <- countTable[,colnames(countTable) %in% row.names(coldata_DGF)]
countTable_DGF <- countTable_DGF[,match(row.names(coldata_DGF), colnames(countTable_DGF))]
model <- as.formula( ~ DSEX + Sample_Set + PC1Scaled + DGF_bin)
dds <- DESeqDataSetFromMatrix(countData = countTable_DGF,
                              colData = coldata_DGF,
                              design = model)
dds <- dds [rowSums(counts(dds))>0,]
dds <- DESeq(dds, betaPrior = FALSE)
saveRDS(dds, file = "DGF_dds.rds")
res <- results(dds, alpha = 0.05)
res <- as.data.frame(res) %>% 
  filter(padj <0.05)
# 2. EGFR3
coldata_EGFR3 <- coldata %>%
  filter(!is.na(EGFR3)) %>% 
  mutate(EGFR3_cont = scale(EGFR3))
countTable_EGFR3 <- countTable[,colnames(countTable) %in% row.names(coldata_EGFR3)]
countTable_EGFR3 <- countTable_EGFR3[,match(row.names(coldata_EGFR3), colnames(countTable_EGFR3))]
model <- as.formula( ~ DSEX + Sample_Set + PC1Scaled + EGFR3_cont)
dds <- DESeqDataSetFromMatrix(countData = countTable_EGFR3,
                              colData = coldata_EGFR3,
                              design = model)
dds <- dds [rowSums(counts(dds))>0,]
dds <- DESeq(dds, betaPrior = FALSE)
saveRDS(dds, file = "EGFR3_dds.rds")
res <- results(dds, alpha = 0.05)
res <- as.data.frame(res) %>% 
  filter(padj <0.05)
# 3. EGFR12
coldata_EGFR12 <- coldata %>%
  filter(!is.na(EGFR12)) %>% 
  mutate(EGFR12_cont = scale(EGFR12))
countTable_EGFR12 <- countTable[,colnames(countTable) %in% row.names(coldata_EGFR12)]
countTable_EGFR12 <- countTable_EGFR12[,match(row.names(coldata_EGFR12), colnames(countTable_EGFR12))]
model <- as.formula( ~ DSEX + Sample_Set + PC1Scaled + EGFR12_cont)
dds <- DESeqDataSetFromMatrix(countData = countTable_EGFR12,
                              colData = coldata_EGFR12,
                              design = model)
dds <- dds [rowSums(counts(dds))>0,]
dds <- DESeq(dds, betaPrior = FALSE)
saveRDS(dds, file = "EGFR12_dds.rds")
res <- results(dds, alpha = 0.05)
res <- as.data.frame(res) %>% 
  filter(padj <0.05)
# 4. Rejection
coldata_rej <- coldata %>%
  mutate(rejection = ifelse(REJECTION_COUNT12 >0, "rejection",
                        ifelse(REJECTION_COUNT3>0, "rejection",
                               ifelse(REJECTION_COUNT12==0, "no_rejection", NA))))%>% 
  filter(!is.na(rejection))
countTable_rej <- countTable[,colnames(countTable) %in% row.names(coldata_rej)]
countTable_rej <- countTable_rej[,match(row.names(coldata_rej), colnames(countTable_rej))]
model <- as.formula( ~ DSEX + Sample_Set + PC1Scaled + rejection)
dds <- DESeqDataSetFromMatrix(countData = countTable_rej,
                              colData = coldata_rej,
                              design = model)
dds <- dds [rowSums(counts(dds))>0,]
dds <- DESeq(dds, betaPrior = FALSE)
saveRDS(dds, file = "rej_dds.rds")
res <- results(dds, alpha = 0.05)
res <- as.data.frame(res) %>% 
  filter(padj <0.05)
