#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

library(WGCNA) # V1.71
library(cluster) # V2.1.4
library(flashClust) # V1.01-2
library(tidyverse) # V1.3.2
library(DESeq2) # V1.34.0
library(limma) # V3.50.3

set.seed(1255342)
options(stringsAsFactors = F)

# load countTable and coldata
coldata <- read.csv("Coldata.csv")
row.names(coldata) <- coldata$X
coldata$X <- NULL
countTable <- read.csv("CountTable.csv")
row.names(countTable) <- countTable$X
countTable$X <- NULL

# Filter for highly cortical samples only
coldata <- filter(coldata, Tissue == "cortex")
countTable <- countTable[, match(row.names(coldata), colnames(countTable))]

coldata[,'Sample_Set']<-factor(coldata[,'Sample_Set'])
coldata[,'Cluster']<-factor(coldata[,'Cluster'])
coldata <- droplevels(coldata)

## remove NAs
coldata <- coldata %>% 
  filter(!is.na(DGF_days)) %>%
  filter(!is.na(EGFR3)) %>%  
  filter(!is.na(EGFR12))
countTable <- countTable[, match(row.names(coldata), colnames(countTable))]

# create outcome variables
coldata <- coldata %>% 
  mutate(DGF_bin = cut(DGF_days, breaks = c(-0.5, 1.5, Inf), 
                       labels = c("noDGF", "DGF"))) %>% 
  mutate(EGFR3_cont = scale(EGFR3)) %>% 
  mutate(EGFR12_cont = scale(EGFR12)) %>% 
  mutate(DAGE_cont = scale(DAGE)) %>% 
  mutate(DBMI_cont = scale(DBMI))

## binarise variables
tmp1= binarizeCategoricalVariable(coldata$DSEX, includePairwise = TRUE, includeLevelVsAll = FALSE)
tmp2= binarizeCategoricalVariable(coldata$DGF_bin, includePairwise = TRUE, includeLevelVsAll = FALSE)
coldata <- data.frame(coldata, tmp1, tmp2)
rm(tmp1, tmp2)

#prepare countTable for WGCNA
model <- as.formula( ~ Sample_Set + DSEX + PC1Scaled + DGF_bin + EGFR12_cont)
dds <- DESeqDataSetFromMatrix(countData = countTable,
                              colData = coldata,
                              design = model)
dds <- DESeq(dds, betaPrior = FALSE)
## Prefilter lowly expressed genes
dds <- dds [ rowSums(counts(dds)) > 50,]
## apply variance stabilising transformation
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)
# Remove batch effects
mm <- model.matrix(~ MALE.vs.FEMALE + PC1Scaled + DGF_bin + EGFR12_cont, colData(vsd))
mat <- limma::removeBatchEffect(mat, batch=vsd$Sample_Set, design=mm)
# Transpose
mat <- t(mat)
# Match row names
coldata <- coldata[match(row.names(mat), row.names(coldata)),]
rm(vsd, mm, model, countTable)


### samples analysis and removal of outliers
## sample adjacency matrix (A), connectivity (k) and standardised connectivity (Z.k)
A.sample <- adjacency(t(mat), type = "distance")
k.sample <- as.numeric(apply(A.sample, 2, sum)) - 1
Z.k.sample = scale(k.sample)
## set Z.k threshold for outlier
threshold_outlier <- -5
outlierC <- ifelse(Z.k.sample < threshold_outlier, "red", "black")
## cluster samples
sampleTree <- flashClust(as.dist(1 - A.sample), method = "average")
## Convert traits to colour representation where red = high values
## select variables of interest
## select variables of interest
coldata <- coldata %>% 
  select(c("DAGE_cont", 
           "DBMI_cont", 
           "MALE.vs.FEMALE",
           "PC1Scaled",
           "DGF.vs.noDGF",
           "EGFR3_cont", 
           "EGFR12_cont"))
traitColors <- data.frame(numbers2colors(coldata, signed = FALSE))
dimnames(traitColors)[[2]] <- paste(names(coldata), "C", sep = "")
datColors <- data.frame(outlierC, traitColors)
pdf(file="sample_tree.pdf", width = 14, height = 7)
plotDendroAndColors(sampleTree, 
                    groupLabels = names(datColors),
                    colors = datColors,
                    main = "Sample dendrogram and trait heatmap") 
dev.off()
## if there are any outliers, remove them
if(sum(outlierC != "black")>0){
  # Remove outlying samples from expression and trait data
  remove.samples <-  Z.k.sample < thresholdZ.k | is.na(Z.k.sample)
  mat <- mat[-remove.samples,]
  coldata <-  coldata[-remove.samples,]
  # Recompute the sample network among the remaining samples
  A.sample <- adjacency(t(mat),type = "distance")
  # Recompute the Z.k values of outlyingness
  k.sample <- as.numeric(apply(A.sample,2,sum))-1
  Z.k.sample <- scale(k.sample)
  rm(remove.samples)
}
## clear sample objects
rm(sampleTree, A.sample, datColors, outlierC, traitColors, Z.k.sample, k.sample, threshold_outlier)
saveRDS(mat, file = "mat.rds")

### Choose power threshold
powers <- c(1:20) 
sft <- pickSoftThreshold(mat, powerVector = powers, networkType = "signed")
pdf(file="threshold_plots.pdf", width = 14, height = 7)
par(mfrow=c(1,2))
# SFT index as a function of different powers
plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n",main=paste("Scale independence"))
text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1],sft$fitIndices[,5],type="n",
     xlab="Soft Threshold (power)",ylab="Mean Connectivity",main=paste("Mean connectivity"))
text(sft$fitIndices[,1],sft$fitIndices[,5],labels=powers,col="red")
dev.off()
# choose power 12 to satisfy R^2 >0.9
power = sft$fitIndices %>% 
  filter(SFT.R.sq >0.9) %>% 
  pull(Power) %>% 
  min()


### Create network
## calculate weighted adjacency matrix - signed, bicor
A <-  adjacency(mat, power = power, type="signed", corFnc = "bicor") 
saveRDS(A, file = "A.rds")
## define a dissimilarity based on the topological overlap
dissTOM <-  TOMdist(A, TOMType = "signed") 
saveRDS(dissTOM, file = "dissTOM.rds")
### Create clusters
## create hierarchical clustering tree
geneTree <-  flashClust(as.dist(dissTOM), method="average")
saveRDS(geneTree, file = "geneTree.rds")
moduleLabels1 <- cutreeDynamic(dendro = geneTree,
                               distM = dissTOM,
                               method = "hybrid",
                               deepSplit = 0,
                               pamRespectsDendro = F,
                               minClusterSize = 30
)
moduleColors1 <- labels2colors(moduleLabels1) 


### Find eigengenes
MEList1 <- moduleEigengenes(mat, colors = moduleColors1) 
MEs1 <-  MEList1$eigengenes
## plot eigengene network/correlations
pdf(file="eigengenes_premerge.pdf", width = 14, height = 7)
plotEigengeneNetworks(MEs1,
                      "",
                      marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),
                      cex.lab=0.8,
                      xLabelsAngle=90) 
dev.off()
### Merge eigengenes (automatic) which are highly correlated
mergingThresh <- 0.20 # default
merged <- mergeCloseModules(mat, moduleColors1, cutHeight = mergingThresh)
saveRDS(merged, file = "merged.rds")
## resulting merged module colors
moduleColors2 <-  merged$colors
saveRDS(moduleColors2, file = "newmodulecolors.rds")
## eigengenes of the newly merged modules:
MEs2 <-  merged$newMEs
# Plot the relationships among the new merged eigengenes
pdf(file="eigengenes_merged.pdf", width = 14, height = 7)
plotEigengeneNetworks(MEs2,
                      "",
                      marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),
                      cex.lab=0.8,
                      xLabelsAngle=90) 
dev.off()

### plot dendrogram and colours
## create colours for GS for variables of interest
# variables of interest: "DAGE", "DBMI", "DGF", "EGFR3_cont", "EGFR12_cont", "PC1Scaled"
variables_interest <- colnames(coldata)
varColors <- data.frame(matrix(vector(),ncol(mat),0),stringsAsFactors = FALSE)
for(i in 1:(length(variables_interest))){
  var.name <- as.character(variables_interest[i])
  var <- as.data.frame(select(coldata, var.name))
  names(var) <-  as.character(variables_interest[i])
  # Next use this trait to define a gene significance variable
  tmp <- as.numeric(cor(mat, var, use="p"))
  # This translates the numeric values into colors
  tmp <- numbers2colors(tmp,signed = T)
  varColors <- cbind(varColors, tmp)
}
## combine all colour data
datColors <- data.frame(moduleColors1, moduleColors2, varColors)
## plot dendrogram (with above colours)
pdf(file="dendrogram.pdf", width = 14, height = 7)
plotDendroAndColors(geneTree,
                    colors = datColors,
                    groupLabels=c("modules1", "modulesmerged", colnames(coldata)),
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05) 
dev.off()
pdf(file="dendrogram2.pdf", width = 14, height = 7)
plotDendroAndColors(geneTree,
                    colors = datColors[,2],
                    groupLabels=c("modulesmerged"),
                    dendroLabels=FALSE,
                    hang=0.03,
                    addGuide=TRUE,
                    guideHang=0.05) 
dev.off()


### link modules with physiological traits (Langfelder tombstone plot)
## Choose a module assignment
mColors <- moduleColors2
## Define numbers of genes and samples
nGenes <-  ncol(mat)
nSamples <-  nrow(mat)
## Recalculate MEs with color labels
MEs <-  moduleEigengenes(mat,mColors)$eigengenes
MEsOrdered = orderMEs(MEs)
modTraitCor = cor(MEsOrdered, coldata, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)
#Since we have a moderately large number of modules and traits,
#a suitable graphical representation will help in reading
#the table. We color code each association by the correlation value:
# Will display correlations and their p-values
textMatrix = paste(signif(modTraitCor, 2), "\n(",
                   signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
par(mar = c(6, 8.5, 3, 3))
# Display the correlation values within a heatmap plot
pdf(file="tombstone_plot.pdf", width = 12, height = 8)
labeledHeatmap(Matrix = modTraitCor, 
               xLabels = names(coldata),
               yLabels = names(MEsOrdered), 
               ySymbols = names(MEsOrdered),
               colorLabels =FALSE,
               colors=blueWhiteRed(50),
               textMatrix=textMatrix,
               setStdMargins = FALSE, 
               cex.text = 0.75, 
               zlim = c(-1,1),
               main = paste("Module-trait relationships")
) 
dev.off()

# gene colours
genecolors <- as.data.frame(colnames(mat))
genecolors$colour <- moduleColors2
colnames(genecolors)[1] <- "gene"
saveRDS(genecolors, file = "genecolors.rds")

### Gene significance measures
## define variables of interest
DGF <- as.data.frame(coldata$DGF.vs.noDGF)
names(DGF) <- "DGF"
EGFR3 <- as.data.frame(coldata$EGFR3_cont)
names(EGFR3) <- "EGFR3"
EGFR12 <- as.data.frame(coldata$EGFR12_cont)
names(EGFR12) <- "EGFR12"
DAGE <- as.data.frame(coldata$DAGE_cont)
names(DAGE) <- "DAGE"
PC1Scaled <- as.data.frame(coldata$PC1Scaled)
names(PC1Scaled) <- "PC1Scaled"
## define gene significance
GS.DGF <- as.numeric(cor(mat, DGF, use = "p"))
GS.EGFR3 <- as.numeric(cor(mat, EGFR3, use = "p"))
GS.EGFR12 <- as.numeric(cor(mat, EGFR12, use = "p"))
GS.DAGE <- as.numeric(cor(mat, DAGE, use = "p"))
GS.PC1Scaled <- as.numeric(cor(mat, PC1Scaled, use = "p"))
# save gene significance
GS.scores <- matrix(ncol = 3, nrow = length(GS.DGF))
row.names(GS.scores) <- colnames(mat)
GS.scores[,1] <- GS.DGF
GS.scores[,2] <- GS.EGFR3
GS.scores[,3] <- GS.EGFR12
colnames(GS.scores) <- c("DGF", "EGFR3", "EGFR12")
saveRDS(GS.scores, file = "GS.scores.rds")
## translate numerics into colours
GF.DGF.colour <- numbers2colors(GS.DGF, signed = T) 
GS.EGFR3.colour <- numbers2colors(GS.EGFR3, signed = T)
GS.EGFR12.colour <- numbers2colors(GS.EGFR12, signed = T)
GS.DAGE.colour <- numbers2colors(GS.DAGE, signed = T)
GS.PC1Scaled.colour <- numbers2colors(GS.PC1Scaled, signed = T)


### Module membership using kME
datKME <- signedKME(mat, MEs2)
saveRDS(datKME, file = "datKME.rds")

### Intramodular analysis
selectModules <- c("orange",
                   "lightcyan") # 
## define interesting outcomes
selectOutcomes <- c("DAGE",
                    "DGF",
                    "EGFR3",
                    "EGFR12")
colorOfColumn <- substring(names(datKME),4)
par(mfrow = c(2,2))
par(mfrow=c(2,length(selectModules)/2))
# plot GS vs. kME for all modules/outcomes
for (i in 1:length(selectOutcomes)){
  tmp <- selectOutcomes[i]
  pdf(file=paste("GSvskME_", tmp, ".pdf", sep = ""), width = 14, height = 7)
  for (module in selectModules) {
    column = match(module,colorOfColumn)
    restModule = mColors == module
    verboseScatterplot(datKME[restModule,column],get(paste("GS.", tmp, sep = ""))[restModule],
                       xlab=paste("Module Membership ",module,"module"),ylab=paste("GS.", tmp, sep = ""),
                       main=paste("kME.",module,"vs. GS"),col=module)
  }
  dev.off()
}

### GSEA for interesting modules
## rank genes by kME values for interesting modules
## create .rnk files for all interesting modules
for (i in 1:length(selectModules)){
  module <- selectModules[i]
  tmp1 <- datKME %>% 
    dplyr::select(paste("kME", module, sep = "")) %>% 
    tibble::rownames_to_column(var = "gene")
  tmp1 <- mutate(tmp1, GSEAmetric = tmp1[,2])
  tmp1 <- tmp1[, c(1, 3)]
  # check for duplicates
  if (sum(duplicated(tmp1$gene))!=0) {
    tmp1 <- tmp1[complete.cases(tmp1),]
  }
  # export as tsv
  filename <- paste("GSEA_", module, ".rnk", sep = "")
  write.table(tmp1, file = filename, quote = F, sep = "\t", row.names = F)
}

### Run all GSEA through GSEA software (V4.1.0)
# Settings: No_collapse, classic, default settings

## Hallmarks
MasterHallmarks <- data.frame(matrix(vector(),0,4,dimnames = list(c(), c("NES", "FDR.q.val", "Comparision", "pathway"))),stringsAsFactors = FALSE)
## specify files
HF <- list.files("Hallmarks/")
## build GSEA plots
for (x in 1:length(HF)){
  tmp <- list.files(paste("Hallmarks/",HF[x],"/", sep=""))
  tmp <- grep("gsea_report", tmp, value = TRUE)
  tmp <- grep("tsv", tmp, value = TRUE)
  tmp <- paste("Hallmarks/",HF[x],"/",tmp, sep="")
  tmp1 <- read.table(tmp[1], fill = TRUE, sep="\t", header=TRUE)
  tmp2 <- read.table(tmp[2], fill = TRUE, sep="\t", header=TRUE)
  Hallmarks <- rbind(tmp1, tmp2)
  tmp <- substr(Hallmarks$NAME, 10,100)
  tmp <- gsub("_", " ", tmp)
  row.names(Hallmarks) <- tmp
  Hallmarks <- Hallmarks[,c(5,6,8)]
  tmp <- strsplit(HF[x],"Gsea")
  tmp <-  as.character(lapply(tmp, function(x) {x[1]}))
  tmp <-substr(tmp, 1, nchar(tmp)-1)
  assign(tmp, Hallmarks)
  Hallmarks$Comparision <- tmp
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
            coord_flip()+
            ggtitle(tmp))
    dev.off()
  }  else if (all(Hallmarks$NES < 0)){
    ## If all hallmarks are negative
    pdf(file=paste("Hallmarks_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
    print(ggplot( data=Hallmarks, aes(x=Hallmark, y=NES, size=-log10(FDR.q.val+0.00001))) +
            geom_point(colour="dodgerblue2") +
            geom_hline(yintercept=0, colour="grey50") +
            labs(y="NES", x="Hallmarks", size="-log10(FDR+0.0001)", colour="")+
            theme_bw()+
            coord_flip()+
            ggtitle(tmp)
    )
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
            coord_flip()+
            ggtitle(tmp))
    dev.off()
  }
}

## customsig plots
MasterCustomsigs <- data.frame(matrix(vector(),0,4,dimnames = list(c(), c("NES", "FDR.q.val", "Comparision", "pathway"))),stringsAsFactors = FALSE)
## specify files
HF <- list.files("Customsigs/")
## build GSEA plots
for (x in 1:length(HF)){
  tmp <- list.files(paste("Customsigs/",HF[x],"/", sep=""))
  tmp <- grep("gsea_report", tmp, value = TRUE)
  tmp <- grep("tsv", tmp, value = TRUE)
  tmp <- paste("Customsigs/",HF[x],"/",tmp, sep="")
  tmp1 <- read.table(tmp[1], fill = TRUE, sep="\t", header=TRUE)
  tmp2 <- read.table(tmp[2], fill = TRUE, sep="\t", header=TRUE)
  Customsigs <- rbind(tmp1, tmp2)
  tmp <- Customsigs$NAME
  tmp <- gsub("_", " ", tmp)
  row.names(Customsigs) <- tmp
  Customsigs <- Customsigs[,c(5,6,8)]
  tmp <- strsplit(HF[x],"Gsea")
  tmp <-  as.character(lapply(tmp, function(x) {x[1]}))
  tmp <-substr(tmp, 1, nchar(tmp)-1)
  assign(tmp, Customsigs)
  Customsigs$Comparision <- tmp
  Customsigs$cell <- row.names(Customsigs)
  MasterCustomsigs <- rbind(MasterCustomsigs, Customsigs)
  ## Remove any sets which are non signifant (FDR > 0.05)
  Customsigs <- Customsigs [Customsigs$FDR.q.val< 0.05,]
  Customsigs <- Customsigs[order(Customsigs$NES, decreasing = FALSE),]
  Customsigs$Hallmark <- factor(row.names(Customsigs), levels = unique(row.names(Customsigs)))
  Customsigs <- droplevels(Customsigs)
  Customsigs$dir <- factor(Customsigs$NES > 0, levels = c("FALSE", "TRUE"))
  if (all(Customsigs$NES > 0)) {
    ## If all Hallmarks are positve
    pdf(file=paste("Customsigs_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
    print(ggplot( data=Customsigs, aes(x=Hallmark, y=NES,  size=-log10(FDR.q.val+0.00001))) +
            geom_point(colour="red") +
            geom_hline(yintercept=0, colour="grey50") +
            labs(y="NES", x="Customsigs", size="-log10(FDR+0.0001)", colour="")+
            theme_bw()+
            coord_flip()+
            ggtitle(tmp))
    dev.off()
  }  else if (all(Customsigs$NES < 0)){
    ## If all hallmarks are negative
    pdf(file=paste("Customsigs_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
    print(ggplot( data=Customsigs, aes(x=Hallmark, y=NES, size=-log10(FDR.q.val+0.00001))) +
            geom_point(colour="dodgerblue2") +
            geom_hline(yintercept=0, colour="grey50") +
            labs(y="NES", x="Customsigs", size="-log10(FDR+0.0001)", colour="")+
            theme_bw()+
            coord_flip()+
            ggtitle(tmp)
    )
    dev.off()
  }  else {
    ## Plot bar graphs
    pdf(file=paste("Customsigs_Alt_",tmp,".pdf",sep = ""), useDingbats = FALSE)
    print(ggplot( data=Customsigs, aes(x=Hallmark, y=NES, colour=dir, size=-log10(FDR.q.val+0.00001))) +
            geom_point() +
            geom_hline(yintercept=0, colour="grey50") +
            labs(y="NES", x="Customsigs", size="-log10(FDR+0.0001)", colour="")+
            scale_colour_manual(breaks = c("FALSE","TRUE"), values=c("dodgerblue2","red"), labels=c("Downregulated\nPathway", "Upregulated\nPathway"))+
            theme_bw()+
            coord_flip()+
            ggtitle(tmp))
    dev.off()
  }
}

# Enrichment plots and leading edge analysis
library(cowplot) # V1.1.1
file_locations <- c("Hallmarks/lightcyan.GseaPreranked.1651416210113/",
                    "Hallmarks/orange.GseaPreranked.1651416201003/")
pathways_lightcyan <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MITOTIC_SPINDLE"
)
pathways_orange <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_MITOTIC_SPINDLE"
)
pathways_all <- list(pathways_lightcyan,
                     pathways_orange)
analyses <- c("lightcyan", 
              "orange")

## Make enrichment plots and find LE genes
for(i in 1: length(pathways_all)){
  pathways <- pathways_all[[i]]
  file_location <- file_locations[i]
  analysis <- analyses[i]
  for (x in 1:length(pathways)){
    tmp1 <- read.table(file = paste(file_location, pathways[x], ".tsv", sep = ""), sep = '\t', header = T)
    p1 <- ggplot(data=tmp1, aes(x=RANK.IN.GENE.LIST, y=RUNNING.ES))+
      geom_line()+
      geom_hline(yintercept=0, color='black')+
      geom_vline(xintercept = 0, color='black') +
      xlab("Gene Rank")+
      ylab("Running Enrichment Score")+
      theme_gray()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())
    ## Remove 1st and last row
    tmp2 <- tmp1[c(-1, -nrow(tmp1)),]
    col<- if (mean(tmp2$RUNNING.ES)> 0)"red" else "dodgerblue2"
    p2 <- ggplot(data=tmp2, aes(y=RANK.IN.GENE.LIST, x=1))+
      geom_violin(fill=col)+
      theme(axis.text.x = element_blank())+
      ylab("Gene Rank")+
      geom_hline(yintercept = 0)+
      theme_gray()+
      expand_limits(y=c(0,nrow(tmp1)))+
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
      coord_flip()
    pdf(file = paste("EP_", analysis, "_", pathways[x], ".pdf", sep = ""), useDingbats = FALSE)
    print(plot_grid(p1, p2, ncol=1, rel_heights = c(4,1), align = "v"))
    dev.off()
    tmp3 <- tmp1 %>% 
      dplyr::filter(CORE.ENRICHMENT == "Yes") %>% 
      dplyr::arrange(desc(abs(RANK.METRIC.SCORE)))
    # write leading edge genes
    write.table(tmp3, file = paste("LE_", analysis, "_", pathways[x], ".tsv", sep = ""), quote=FALSE, sep='\t', col.names = T)
  }
}
