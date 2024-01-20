#  R version 4.1.3 (2022-03-10) -- "One Push-Up"

library(tidyverse) # V1.3.2
library(DESeq2) # V1.34.0
library(ggrepel) # V0.9.2
library(pheatmap) # V1.0.12
library(cowplot) # V1.1.1

set.seed(1255342)

# pathway analysis

# Import DE results
DGF_dds <- readRDS(file = "DGF_dds.rds")
EGFR3_dds <- readRDS(file = "EGFR3_dds.rds")
EGFR12_dds <- readRDS(file = "EGFR12_dds.rds")
rej_dds <- readRDS(file = "rej_dds.rds")
dds_all <- list(DGF_dds, EGFR3_dds, EGFR12_dds, rej_dds)
results <- function(x){
  DESeq2::results(x, alpha = 0.05)
}
res_all <- lapply(dds_all, results)
analyses <- c("DGF_bin", 
              "EGFR3_cont", 
              "EGFR12_cont",
              "rej")

# GSEA
# create pre-ranked files
for(i in 1:(length(res_all))){
  res <- res_all[[i]]
  # Creating ranked gene list 
  GSEArnk <- as.data.frame(na.omit(res))
  GSEArnk <- rownames_to_column(GSEArnk, var = "gene")
  # rank by 1/pvalue*sign(logfold)
  GSEArnk <- GSEArnk %>% 
    mutate(GSEAmetric = 1/GSEArnk$pvalue * sign(GSEArnk$log2FoldChange)) %>% 
    arrange(desc(GSEAmetric)) %>% 
    dplyr::select(c("gene", "GSEAmetric"))
  # check for duplicates
  if (sum(duplicated(GSEArnk$gene))!=0) {
    GSEArnk <- GSEArnk[complete.cases(GSEArnk),]
  }
  # export as tsv
  filename <- paste("GSEA_", analyses[i], ".rnk", sep = "")
  write.table(GSEArnk, file = filename, quote = F, sep = "\t", row.names = F)
}

### Run all GSEA through GSEA software (V4.1.0)
# Settings: No_collapse, classic, default settings

# Hallmarks
MasterHallmarks <- data.frame(matrix(vector(),0,4,dimnames = list(c(), c("NES", "FDR.q.val", "Comparision", "pathway"))),stringsAsFactors = FALSE)
HF <- list.files("Hallmarks/") # folder for output files

# build GSEA plots
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
file_locations <- c("Hallmarks/DGF.GseaPreranked.1649467694105/",
                    "Hallmarks/EGFR3.GseaPreranked.1649467726739/",
                    "Hallmarks/EGFR12.GseaPreranked.1649467760123/")
pathways_DGF <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_MITOTIC_SPINDLE"
)
pathways_EGFR3 <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)
pathways_EGFR12 <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MYOGENESIS",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_PROTEIN_SECRETION",
  "HALLMARK_MTORC1_SIGNALING"
)
pathways_all <- list(pathways_DGF,
                     pathways_EGFR3,
                     pathways_EGFR12)
analyses <- c("DGF_bin", 
              "EGFR3_cont", 
              "EGFR12_cont")

for(i in 1: length(analyses)){
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

# plot heatmaps for top n leading edge genes in each pathway of interest
# file names for LE genes
tmp <- list.files(".")
LEs_DGF <- grep("LE_DGF_bin", tmp, value = TRUE)
LEs_EGFR3 <- grep("LE_EGFR3_cont", tmp, value = TRUE)
LEs_EGFR12 <- grep("LE_EGFR12_cont", tmp, value = TRUE)
LEs_all <- list(LEs_DGF, LEs_EGFR3, LEs_EGFR12)
dds_all <- list(DGF_dds, EGFR3_dds, EGFR12_dds)
# specify top n leading edge genes
number = 20
# plot heatmaps
for(x in 1:length(LEs_all)){
  LEs <- LEs_all[[x]]
  analysis <- analyses[x]
  dds <- dds_all[[x]]
  for(i in 1:length(LEs)){
    filename <- LEs[i]
    tmp <- read.table(file = paste(filename), sep = "\t")
    tmp <- tmp %>% 
      slice_head(n = number)
    tmp2 <- tmp$SYMBOL
    ## Extract transformed counts for these genes
    plotDat <- vst(dds, blind = FALSE)[tmp2,] %>% 
      assay()
    # Set annotations
    anno_info <- as.data.frame(colData(dds)[,c("Sample", analysis)])
    # Change order of samples
    anno_info <- anno_info[order(anno_info[, 2]),]
    plotDat <- plotDat[,match(row.names(anno_info), colnames(plotDat))]
    anno_info$Sample <- NULL
    summary(colnames(plotDat)==row.names(anno_info))
    tmp3 <- sapply(strsplit(filename,"_HALLMARK_"), `[`, 2)
    tmp3 <- substr(tmp3,1,nchar(tmp3) - 4)
    filename <- paste("LEheatmap_", analysis, "_HALLMARK_", tmp3, ".pdf", sep = "")
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
    pdf(file=filename, useDingbats = FALSE)
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
    dev.off()
  }
}


# Create gene lists for String db
# Top n up and down
number <- 50
for (i in 1:(length(res_all))){
  res <- res_all[[i]]
  res <- as.data.frame(res)
  res_up <- res %>% 
    rownames_to_column(var = "gene") %>% 
    filter(log2FoldChange > 0) %>% 
    arrange(padj) %>% 
    slice_head(n = number)
  filename <- paste("stringdb_", analyses[i], "_top", number, "up.tsv", sep = "")
  write.table(res_up$gene, file = filename, quote = FALSE, sep = '\t', row.names = F, col.names = F)
  res_down <- res %>% 
    rownames_to_column(var = "gene") %>% 
    filter(log2FoldChange < 0) %>% 
    arrange(padj) %>% 
    slice_head(n = number)
  filename <- paste("stringdb_", analyses[i], "_top", number, "down.tsv", sep = "")
  write.table(res_down$gene, file = filename, quote = FALSE, sep = '\t', row.names = F, col.names = F)
}

# visit stringdb.org

