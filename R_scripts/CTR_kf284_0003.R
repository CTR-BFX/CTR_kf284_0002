#!/usr/local/bin/Rscript

#------------------------------------------------------------------------------
# RNA-Seq Analysis to accompany:
#
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/
#
# CTR Code: CTR_kf284_0003
#
# Analysis Performed by Russell S. Hamilton
# CTR Bioinformatics 
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------


#
# initial install of packages
#
#source("http://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
#devtools::install_github('slowkow/ggrepel')
#install.packages("ggrepel")
#install.packages("gplots")

library('DESeq2')
library('ggplot2')
library('RColorBrewer')
library("cowplot")
library("pheatmap")
library("ggrepel")
library("reshape2")
library("biomaRt")
library("matrixStats")
library("plyr")
library("dplyr")
library("UpSetR")
library("genefilter")
library("scales")
library("gtools")
library("Cairo")
library("GenomicRanges")
library("rtracklayer")
library("BiocParallel")

register(MulticoreParam(2))


message("+-------------------------------------------------------------------------------")
message("+ Set up some constants e.g. base directories")
message("+-------------------------------------------------------------------------------")

Project         <- "CTR_kf284_0003"
#Base.dir        <- getwd() 
Base.dir        <- "/Users/rhamilto/Documents/CTR-Groups/Kristian_Franze/CTR_kf284_0003_Mouse"
setwd(Base.dir)
HTSeq.dir       <- paste(Base.dir,"/HTSEQ-COUNTS", sep="")

elementTextSize <- 10

significance    <- 0.01
log2foldchange  <- 2
topN            <- 30

message("+-------------------------------------------------------------------------------")
message("+ Set up the sample table")
message("+-------------------------------------------------------------------------------")

sampleFiles      <- grep('*htseq_counts.txt',list.files(HTSeq.dir),value=TRUE)
sampleBarcodes   <- gsub(".H23LKBBXY.s_3.r_1_trimmed_GRCm38.star.bam_htseq_counts.txt", "", sampleFiles)
sampleBarcodes   <- gsub("SLX-[0-9]*.", "", sampleBarcodes)

sampleCondition  <- sampleBarcodes
sampleCondition  <- gsub("D701rna_D501rna", "PLAIN", sampleCondition)
sampleCondition  <- gsub("D702rna_D501rna", "PLAIN", sampleCondition)
sampleCondition  <- gsub("D703rna_D501rna", "PLAIN", sampleCondition)
sampleCondition  <- gsub("D704rna_D501rna", "PLAIN", sampleCondition)
sampleCondition  <- gsub("D705rna_D501rna", "PLAIN", sampleCondition)
sampleCondition  <- gsub("D706rna_D501rna", "PLAIN", sampleCondition)
sampleCondition  <- gsub("D701rna_D502rna", "DEX",   sampleCondition)
sampleCondition  <- gsub("D702rna_D502rna", "DEX",   sampleCondition)
sampleCondition  <- gsub("D703rna_D502rna", "DEX",   sampleCondition)
sampleCondition  <- gsub("D704rna_D502rna", "DEX",   sampleCondition)
sampleCondition  <- gsub("D705rna_D502rna", "DEX",   sampleCondition)
sampleCondition  <- gsub("D706rna_D502rna", "DEX",   sampleCondition)
sampleCondition  <- gsub("D701rna_D503rna", "NLRP3", sampleCondition)
sampleCondition  <- gsub("D702rna_D503rna", "NLRP3", sampleCondition)
sampleCondition  <- gsub("D703rna_D503rna", "NLRP3", sampleCondition)
sampleCondition  <- gsub("D704rna_D503rna", "NLRP3", sampleCondition)
sampleCondition  <- gsub("D705rna_D503rna", "NLRP3", sampleCondition)
sampleCondition  <- gsub("D706rna_D503rna", "NLRP3", sampleCondition)

sampleNames      <- sampleBarcodes
sampleNames      <- gsub("D701rna_D501rna", "01", sampleNames)
sampleNames      <- gsub("D702rna_D501rna", "02", sampleNames)
sampleNames      <- gsub("D703rna_D501rna", "03", sampleNames)
sampleNames      <- gsub("D704rna_D501rna", "04", sampleNames)
sampleNames      <- gsub("D705rna_D501rna", "05", sampleNames)
sampleNames      <- gsub("D706rna_D501rna", "06", sampleNames)
sampleNames      <- gsub("D701rna_D502rna", "07", sampleNames)
sampleNames      <- gsub("D702rna_D502rna", "08", sampleNames)
sampleNames      <- gsub("D703rna_D502rna", "09", sampleNames)
sampleNames      <- gsub("D704rna_D502rna", "10", sampleNames)
sampleNames      <- gsub("D705rna_D502rna", "11", sampleNames)
sampleNames      <- gsub("D706rna_D502rna", "12", sampleNames)
sampleNames      <- gsub("D701rna_D503rna", "13", sampleNames)
sampleNames      <- gsub("D702rna_D503rna", "14", sampleNames)
sampleNames      <- gsub("D703rna_D503rna", "15", sampleNames)
sampleNames      <- gsub("D704rna_D503rna", "16", sampleNames)
sampleNames      <- gsub("D705rna_D503rna", "17", sampleNames)
sampleNames      <- gsub("D706rna_D503rna", "18", sampleNames)

sampleTable      <- data.frame(sampleName=sampleNames, fileName=sampleFiles, 
                               Barcode=sampleBarcodes, condition=sampleCondition)

print(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Removing samples identified by QC")
message("+-------------------------------------------------------------------------------")

sampleTable <- sampleTable[-c(2, 16),]

print(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Retrieve ensEMBL annotations")
message("+-------------------------------------------------------------------------------")
ensembl    <- useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'entrezgene', 'chromosome_name', 'gene_biotype'), mart = ensembl)  
head(ensEMBL2id)
nrow(ensEMBL2id)


ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene', 'description', 'go_id', 'name_1006' ), mart = ensembl)    
head(ensEMBL2id_go)


# Inflammasome_GO_genes_ensembl
# Inflammatory_GO_genes_ensembl
# fibrosis_GO_genes_ensembl
# Fibrosis_myofibroblast
# Negative_Regulation_of_Inflammasomes
# Signalling_Downstream_of_Inflammasomes
# Main
# Gene.Il10

# inflammasome genes:
ensEMBL2id_go_Inflammasome     <- subset(ensEMBL2id_go, grepl("inflammasome", ensEMBL2id_go$name_1006))
ensEMBL2id_go_Inflammasome     <- subset(ensEMBL2id_go_Inflammasome, ensembl_gene_id != "ENSMUSG00000116526")
Inflammasome_GO_terms          <- unique(ensEMBL2id_go_Inflammasome$name_1006)
Inflammasome_GO_genes_external <- unique(ensEMBL2id_go_Inflammasome$external_gene_name)
Inflammasome_GO_genes_ensembl  <- unique(ensEMBL2id_go_Inflammasome$ensembl_gene_id)
length(Inflammasome_GO_genes_ensembl)

# inflammatory genes:
ensEMBL2id_go_Inflammatory     <- subset(ensEMBL2id_go, grepl("inflamm", ensEMBL2id_go$name_1006))
Inflammatory_GO_terms          <- unique(ensEMBL2id_go_Inflammatory$name_1006)
Inflammatory_GO_genes_external <- unique(ensEMBL2id_go_Inflammatory$external_gene_name)
Inflammatory_GO_genes_ensembl  <- unique(ensEMBL2id_go_Inflammatory$ensembl_gene_id)
length(Inflammatory_GO_genes_ensembl)
# Fibrosis genes:
ensEMBL2id_go_Fibrosis         <- subset(ensEMBL2id_go, (ensEMBL2id_go$go_id == "GO:0002248" | 
                                                           ensEMBL2id_go$go_id == "GO:1904596" | 
                                                           ensEMBL2id_go$go_id == "GO:0097709") )
fibrosis_GO_terms              <- unique(ensEMBL2id_go_Fibrosis$name_1006)
fibrosis_GO_genes_external     <- unique(ensEMBL2id_go_Fibrosis$external_gene_name)
fibrosis_GO_genes_ensembl      <- unique(ensEMBL2id_go_Fibrosis$ensembl_gene_id)
length(fibrosis_GO_genes_ensembl)

Fibrosis_myofibroblast <- c("ENSMUSG00000030968","ENSMUSG00000030970","ENSMUSG00000030983","ENSMUSG00000030986","ENSMUSG00000030987",
                            "ENSMUSG00000031012","ENSMUSG00000031021","ENSMUSG00000080657","ENSMUSG00000080673","ENSMUSG00000035783",	
                            "ENSMUSG00000024563","ENSMUSG00000032402","ENSMUSG00000024515","ENSMUSG00000007613","ENSMUSG00000032440",	
                            "ENSMUSG00000002603")

Negative_Regulation_of_Inflammasomes <- c("ENSMUSG00000057329", "ENSMUSG00000007659", "ENSMUSG00000031132", 
                                          "ENSMUSG00000021939", "ENSMUSG00000021270", "ENSMUSG00000020048", 
                                          "ENSMUSG00000022534", "ENSMUSG00000032322", "ENSMUSG00000022024", 
                                          "ENSMUSG00000024401", "ENSMUSG00000022015", "ENSMUSG00000005824",
                                          "ENSMUSG00000026700")

Signalling_Downstream_of_Inflammasomes <- c("ENSMUSG00000055170", "ENSMUSG00000027776", "ENSMUSG00000004296",
                                            "ENSMUSG00000039217", "ENSMUSG00000027398", "ENSMUSG00000024810",
                                            "ENSMUSG00000031392", "ENSMUSG00000018899", "ENSMUSG00000032508",
                                            "ENSMUSG00000029468", "ENSMUSG00000031934", "ENSMUSG00000032487",
                                            "ENSMUSG00000056458", "ENSMUSG00000041135", "ENSMUSG00000032041",
                                            "ENSMUSG00000038393")

Main <- c("ENSMUSG00000027995", "ENSMUSG00000039005", "ENSMUSG00000028163", "ENSMUSG00000030793",
          "ENSMUSG00000069830", "ENSMUSG00000070390", "ENSMUSG00000032691", "ENSMUSG00000037860",
          "ENSMUSG00000025888", "ENSMUSG00000033538", "ENSMUSG00000027398", "ENSMUSG00000039217",
          "ENSMUSG00000026981", "ENSMUSG00000039193")

Gene.Il10 <- "ENSMUSG00000016529"


message("+-------------------------------------------------------------------------------")
message("+ Create ddsHTSeq object")
message("+-------------------------------------------------------------------------------")

ddsHTSeq <- ""
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ condition ) 
colData(ddsHTSeq)


message("+-------------------------------------------------------------------------------")
message("+ Create dds objects")
message("+-------------------------------------------------------------------------------")

dds <- ""
dds <- DESeq(ddsHTSeq)
resultsNames(dds)
colData(dds)
colnames(dds)


message("+-------------------------------------------------------------------------------")
message("+ Create Counts object")
message("+-------------------------------------------------------------------------------")

normCounts                     <- counts(dds,      normalized=TRUE)
normCounts.df                  <- as.data.frame(normCounts)
normCounts.df$ensembl_gene_id  <- rownames(normCounts.df)
normCounts.df.annot            <- merge(normCounts.df,ensEMBL2id, by="ensembl_gene_id" )
head(normCounts.df.annot)
#write.csv(normCounts.df.annot, file=paste(Project, "NormalisedCounts_normCounts.csv", sep=""))


message("+-------------------------------------------------------------------------------")
message("+ Run Transformations")
message("+-------------------------------------------------------------------------------")

rld  <- rlogTransformation(dds, blind=T)

vsd  <- varianceStabilizingTransformation(dds, blind=T)


message("+-------------------------------------------------------------------------------")
message("+ Create PCA Plots")
message("+-------------------------------------------------------------------------------")

customPCA <- function(sampleTBL, RLD, TOPNUM, model) {
  
  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, code=sampleTBL$Barcode, condition=sampleTBL$condition)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, shape=condition, label=scores$sampleName) ) +
    
    geom_point(size = 3, alpha=0.75 ) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    scale_shape_manual(name="condition", values = c(17, 16, 15)) + 
    theme(text = element_text(size=elementTextSize)) 
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=scores$sampleName )) +
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) +
    geom_point(size = 3, alpha=0.75 ) + 
    xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
    scale_color_manual(name="Treatment", values = c("PLAIN"="grey", "DEX"="purple3", "NLRP3"="blue") )  +
    scale_fill_manual(name="Treatment", values=c("PLAIN"="grey", "DEX"="purple3", "NLRP3"="blue")) +
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize)) 
  
  return(list(plt.pca, plt.pca.nl) )
}

pca.plt.rld.50             <- customPCA(sampleTable, assay(rld), 50,  "rld")
pca.plt.rld.100            <- customPCA(sampleTable, assay(rld), 100, "rld")
pca.plt.rld.250            <- customPCA(sampleTable, assay(rld), 250, "rld")
pca.plt.rld.500            <- customPCA(sampleTable, assay(rld), 500, "rld")
pca.plt.rld.1000           <- customPCA(sampleTable, assay(rld), 1000, "rld")
pca.plt.rld.all            <- customPCA(sampleTable, assay(rld), nrow(assay(rld)), "rld")


pca.plt.rld.250


pdf(paste(Project, "_Fig.PCA_cleaned.pdf", sep=""),width=15,height=12)
par(bg=NA)
plot_grid(pca.plt.rld.250[[2]], pca.plt.rld.500[[2]], pca.plt.rld.1000[[2]], pca.plt.rld.all[[2]], nrow=2, ncol=2)
dev.off()

#pdf(paste0(Base.dir, "/", Project, "_Fig.PCA.500.pdf"),width=10,height=7)
#par(bg=NA)
#plot_grid(pca.plt.500[[2]], pca.plt.500[[2]], plt.gender, pca.plt.500[[3]], nrow=2, ncol=2, scale = c(0.75, 1, 1, 0.75))
#dev.off()


message("+-------------------------------------------------------------------------------")
message("+ Perform DE with specific contrasts and Annotate genes and Make VCL & MA plots ")
message("+-------------------------------------------------------------------------------")


functionGetSigDEG <- function(result, Project, Title, significance, foldchange, ensEMBL2id) {
  # Filter by adjusted p-value
  result.df <- as.data.frame(result)
  result.df.sub <- as.data.frame(subset(result, padj <= significance  & abs(log2FoldChange) >= foldchange))  
  # Print out the number of significant hits
  message(paste0("+   ", Title, "        padj    = ", nrow(result.df.sub)))
  # Annotate Genes
  result.df$ensembl_gene_id <- rownames(result.df)
  # Annotate from the live ensEMBL table generated above
  result.ann <- merge(result.df,ensEMBL2id, by="ensembl_gene_id")
  # Tidy up the description field
  result.ann$description            <- gsub("..Source.*", "", result.ann$description)
  # Order the results object
  result.ann <- result.ann[order(result.ann$log2FoldChange,decreasing=TRUE),]
  # Write the results to file
  write.csv(result.ann[order(abs(result.ann$log2FoldChange),decreasing=TRUE),], 
            file=paste(Project, "_DESeq2_DEGs_padj_sig", significance, '_fc', foldchange, '_', Title, '.ann.csv', sep=""))
  
  #write.csv(result.ann[order(abs(result.ann$log2FoldChange),decreasing=TRUE),], 
  #          file=paste(Project, "_DESeq2_DEGs_pval_sig", significance, '_fc', foldchange, '_', Title, '.ann.csv', sep=""))
  return(result.ann)
}


message("+-------------------------------------------------------------------------------")
message("+ Create Volcano Plots")
message("+-------------------------------------------------------------------------------")

functionPlotDEVolcano <- function(results, title, sig_cut, logfc_cut, topN, col_gender, xrange, yrange) {
  
  results       <- as.data.frame(results)
  results$genes <- results$external_gene_name
  
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log(padj), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    #geom_point(alpha=0.75, size=2, colour=ifelse(results$padj<=sig_cut & results$log2FoldChange >= logfc_cut,"red","grey")) +
    
    geom_point(data=subset(results, abs(results$log2FoldChange) < logfc_cut | results$padj > sig_cut), alpha=0.75, size=2, colour="grey") +
    geom_point(data=subset(results, results$padj<=sig_cut & results$log2FoldChange >= logfc_cut),      alpha=0.75, size=2, colour="red") +
    geom_point(data=subset(results, results$padj<=sig_cut & results$log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=2, colour="blue") +
    
    geom_point(data=subset(results, results$chromosome_name == "Y" & results$padj<=sig_cut & col_gender == 1), alpha=0.99, size=2, colour="blue") +
    geom_point(data=subset(results, results$chromosome_name == "X" & results$padj<=sig_cut & col_gender == 1), alpha=0.99, size=2, colour="pink") +
    
    geom_text_repel( data= subset(results, results$log2FoldChange > 0      & results$padj<=sig_cut)[1:topN,], 
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= tail(subset(results, results$log2FoldChange < 0 & results$padj<=sig_cut),topN),    
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab("log2FC") + ylab(bquote("-log"[10]~"(adj.p.value)")) + 
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) + 
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    ggtitle(paste0(title, "\n[adj.pval < ", sig_cut, ", abs(l2fc) >", logfc_cut, ", labelled = top ", topN, "]" )) + theme(aspect.ratio=1)
  
  return(volc.plt)
}

message("+-------------------------------------------------------------------------------")
message("+ Create MA Plots")
message("+-------------------------------------------------------------------------------")

functionCustomMAPlot <- function(results, Project, FigureID, Title, significance, log2FC, topN, col_gender, ylim, yrange) {
  
  results$log2FoldChange[results$log2FoldChange > ylim]    <- ylim
  results$log2FoldChange[results$log2FoldChange < -(ylim)] <- -(ylim)
  
  res_ord_red  <- results[order(results$log2FoldChange,decreasing=TRUE),] 
  res_ord_red  <- subset(res_ord_red, log2FoldChange >= log2FC & padj <= significance)
  
  res_ord_blue <- results[order(results$log2FoldChange,decreasing=FALSE),] 
  res_ord_blue <- subset(res_ord_blue, log2FoldChange <= log2FC & padj <= significance)
  
  plt <- ggplot(data = results, aes(x=baseMean, y=log2FoldChange )) +
    geom_abline(intercept = 0,       slope = 0, colour='black', alpha=0.5) +
    geom_abline(intercept = log2FC,  slope = 0, colour='black', alpha=0.5, linetype="dashed") +
    geom_abline(intercept = -log2FC, slope = 0, colour='black', alpha=0.5, linetype="dashed") +
    geom_point(size=0.5, alpha=0.5, col="black") +
    geom_point(data=subset(results, padj <= significance & log2FoldChange > 0), size=1, alpha=0.95,  col="red") +
    geom_point(data=subset(results, padj <= significance & log2FoldChange < 0),  size=1, alpha=0.95,  col="blue") +
    
    geom_point(data=subset(results, results$chromosome_name == "Y" & padj <= significance & col_gender == 1), alpha=0.99, size=1, colour="lightblue") +
    geom_point(data=subset(results, results$chromosome_name == "X" & padj <= significance & col_gender == 1), alpha=0.99, size=1, colour="pink") +
    
    geom_label_repel(data=res_ord_red[c(1:topN),], aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                     size=3, segment.size = 0.25, segment.color = 'darkred', force=5, nudge_x = 0, nudge_y=0) +
    geom_label_repel(data=res_ord_blue[c(1:topN),],
                     aes( x=baseMean, y=log2FoldChange, label=external_gene_name),
                     fill='white', colour='black', point.padding = unit(0.1, "lines"),  
                     size=3, segment.size = 0.25, segment.color = 'darkblue', force=5, nudge_x = 0, nudge_y=0) +
    scale_x_log10( ) +
    scale_y_continuous(breaks=seq(yrange[1],yrange[2],yrange[3])) + 
    xlab("Mean Normalised Read Count") + ylab("log2 Fold Change") + ggtitle(paste(Project, " ", FigureID, "\n", Title, " [l2fc ", log2FC, ", sig ", significance, "]", sep=""))
  
  return (plt)
}

message("+-------------------------------------------------------------------------------")
message("+ Check and other individual genes")
message("+-------------------------------------------------------------------------------")

makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, TITLE, gene2plot,outdir) {
  #
  # Plot the normalised read counts for a specified gene
  #
  if(missing(outdir)){ outdir = "" }
  else( outdir <- paste(outdir, "/", sep=""))
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")
  
  pdf(paste0(outdir, Project, "-DGE_", gene2plot, "_", TITLE, "_collated.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({
    ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
      geom_violin(trim=TRUE, alpha=.5)  + 
      geom_boxplot(width = 0.1, fill='white') + 
      geom_point(position=position_jitter(w=0.1,h=0), alpha=0.5) +
      ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot)) + 
      xlab("") + ylab("Normalised count") +
      scale_fill_manual(name="Genotype", values = c("purple","purple4", "lightgreen", "green")) +
      theme(text = element_text(size=elementTextSize), legend.position="none") })
  dev.off()
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "condition", "samples")
  t2           <- t2[order(t2$condition),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  pdf(paste0(outdir, Project, "-DGE_", gene2plot, "_", TITLE, "_individual.pdf"),width=10,height=7, onefile=FALSE)
  par(bg=NA)
  print({ ggplot(t2, aes(x=samples2, y=count, fill=condition, group=condition)) + 
      geom_bar(stat="identity", alpha=.5) + 
      scale_fill_manual(name="Comparison", values = c("purple","purple4", "lightgreen", "green")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("") +
      ggtitle(paste0(Project, " ::: ", gene2plot, " ::: ", genename2plot)) })
  dev.off()
  print(paste("Created plot for", gene2plot), sep=" ")
}



resultsNames(dds)

res.condition.DEX.vs.PLAIN   <- lfcShrink(dds=dds, contrast=c("condition","DEX","PLAIN"),   type="normal", parallel=TRUE)
nrow(subset(res.condition.DEX.vs.PLAIN, padj <= significance & abs(log2FoldChange) >= log2foldchange) )
nrow(subset(res.condition.DEX.vs.PLAIN, padj <= significance & log2FoldChange >= log2foldchange) )
nrow(subset(res.condition.DEX.vs.PLAIN, padj <= significance & log2FoldChange <= -(log2foldchange) ))

res.condition.DEX.vs.PLAIN.ann     <- functionGetSigDEG(as.data.frame(res.condition.DEX.vs.PLAIN), Project, "DEX.vs.PLAIN",significance,log2foldchange,ensEMBL2id)
res.condition.DEX.vs.PLAIN.ann.vlc <- functionPlotDEVolcano(res.condition.DEX.vs.PLAIN.ann, "DEX.vs.PLAIN", significance, log2foldchange, 
                                                            topN, col_gender=0, xrange=c(-6,10,2),yrange=c(0,120,20))
res.condition.DEX.vs.PLAIN.ann.ma  <- functionCustomMAPlot(as.data.frame(res.condition.DEX.vs.PLAIN.ann), Project, "MA", "DEX.vs.PLAIN", significance, log2foldchange,
                                                           topN, col_gender=0, ylim=10, yrange=c(-8,10,2))

pdf(paste(Project, "_Fig.DE.res.condition.DEX.vs.PLAIN.pdf", sep=""),width=12,height=5)
par(bg=NA)
plot_grid(res.condition.DEX.vs.PLAIN.ann.vlc, res.condition.DEX.vs.PLAIN.ann.ma, nrow=1, ncol=2, rel_widths=c(0.6,1.0))
dev.off()

makeGeneCountPlot(dds, ensEMBL2id, "condition", "Condition", 'ENSMUSG00000026180') # cxcr2
makeGeneCountPlot(dds, ensEMBL2id, "condition", "Condition", 'ENSMUSG00000026725') # Tnn
makeGeneCountPlot(dds, ensEMBL2id, "condition", "Condition", 'ENSMUSG00000030278') # cidec


res.condition.NLRP3.vs.PLAIN <- lfcShrink(dds=dds, contrast=c("condition","NLRP3","PLAIN"), type="normal", parallel=TRUE)
nrow(subset(res.condition.NLRP3.vs.PLAIN, padj <= significance & abs(log2FoldChange) >= log2foldchange) )
head(subset(res.condition.NLRP3.vs.PLAIN, padj <= significance & abs(log2FoldChange) >= log2foldchange) )

res.condition.NLRP3.vs.PLAIN.ann     <- functionGetSigDEG(as.data.frame(res.condition.NLRP3.vs.PLAIN), Project, "NLRP3.vs.PLAIN",significance,log2foldchange,ensEMBL2id)
res.condition.NLRP3.vs.PLAIN.ann.vlc <- functionPlotDEVolcano(res.condition.NLRP3.vs.PLAIN.ann, "NLRP3.vs.PLAIN", significance, log2foldchange, 
                                                            topN, col_gender=0, xrange=c(-4,4,1),yrange=c(0,20,5))
res.condition.NLRP3.vs.PLAIN.ann.ma  <- functionCustomMAPlot(as.data.frame(res.condition.NLRP3.vs.PLAIN.ann), Project, "MA", "NLRP3.vs.PLAIN", significance, log2foldchange,
                                                           topN, col_gender=0, ylim=10, yrange=c(-8,10,2))

pdf(paste(Project, "_Fig.DE.res.condition.NLRP3.vs.PLAIN.pdf", sep=""),width=12,height=5)
par(bg=NA)
plot_grid(res.condition.NLRP3.vs.PLAIN.ann.vlc, res.condition.NLRP3.vs.PLAIN.ann.ma, nrow=1, ncol=2, rel_widths=c(0.6,1.0))
dev.off()





res.condition.NLRP3.vs.DEX   <- lfcShrink(dds=dds, contrast=c("condition","NLRP3","DEX"),   type="normal", parallel=TRUE)
nrow(subset(res.condition.NLRP3.vs.DEX, padj <= significance & abs(log2FoldChange) >= log2foldchange) )
head(subset(res.condition.NLRP3.vs.DEX, padj <= significance & abs(log2FoldChange) >= log2foldchange) )

res.condition.NLRP3.vs.DEX.ann     <- functionGetSigDEG(as.data.frame(res.condition.NLRP3.vs.DEX), Project, "NLRP3.vs.DEX",significance,log2foldchange,ensEMBL2id)
res.condition.NLRP3.vs.DEX.ann.vlc <- functionPlotDEVolcano(res.condition.NLRP3.vs.DEX.ann, "NLRP3.vs.DEX", significance, log2foldchange, 
                                                              topN, col_gender=0, xrange=c(-8,6,2),yrange=c(0,120,20))
res.condition.NLRP3.vs.DEX.ann.ma  <- functionCustomMAPlot(as.data.frame(res.condition.NLRP3.vs.DEX.ann), Project, "MA", "NLRP3.vs.DEX", significance, log2foldchange,
                                                             topN, col_gender=0, ylim=10, yrange=c(-8,10,2))

pdf(paste(Project, "_Fig.DE.res.condition.NLRP3.vs.DEX.pdf", sep=""),width=12,height=5)
par(bg=NA)
plot_grid(res.condition.NLRP3.vs.DEX.ann.vlc, res.condition.NLRP3.vs.DEX.ann.ma, nrow=1, ncol=2, rel_widths=c(0.6,1.0))
dev.off()




message("+-------------------------------------------------------------------------------")
message("+ KEGG ")
message("+-------------------------------------------------------------------------------")

library("pathview")
library("gage")
library("gageData")
library("dplyr")
library("cowplot")
library("pathview")
library("org.Mm.eg.db")


data(korg)
head(korg[,1:3], n=20)
data(kegg.sets.mm)
data(go.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 1)


customKEGG <- function(RESULTS, OUTFIX, CHOSENKEGG, FULLSET) {
  
  Kegg_genes <- RESULTS[!is.na(RESULTS$entrezgene),]  
  Kegg_genes <- Kegg_genes[order(Kegg_genes$entrezgene, -abs(Kegg_genes$log2FoldChange) ), ]        # sort by id and reverse of abs(value)
  Kegg_genes <- Kegg_genes[!duplicated(Kegg_genes$entrezgene), ]                       # take the first row within each id
  colnames(Kegg_genes)
  colnames(Kegg_genes)[colnames(Kegg_genes)=="entrezgene"] <- "Gene_ID" # exchange column name "entrez" to "Gene_ID"
  
  foldchanges = Kegg_genes$log2FoldChange
  names(foldchanges) = Kegg_genes$Gene_ID
  
  mmu04110 <- pathview(gene.data  = foldchanges,
                       pathway.id = CHOSENKEGG,
                       species    = "mmu", out.suffix = OUTFIX)
  
  
  ego <- enrichGO(gene          = FULLSET,
                  universe      = names(Kegg_genes$Gene_ID),
                  OrgDb         = org.Mm.eg.db,
                  ont           = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.01,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  print(head(ego))
  
  return(ego)
}

customKEGG( subset(res.condition.DEX.vs.PLAIN.ann,   abs(log2FoldChange) >= log2foldchange & padj <= significance), "DEX.vs.PLAIN",   "mmu04621", res.condition.DEX.vs.PLAIN.ann$entrezgene   )
customKEGG( subset(res.condition.NLRP3.vs.PLAIN.ann, abs(log2FoldChange) >= log2foldchange & padj <= significance), "NLRP3.vs.PLAIN", "mmu04621", res.condition.NLRP3.vs.PLAIN.ann$entrezgene )
customKEGG( subset(res.condition.NLRP3.vs.DEX.ann,   abs(log2FoldChange) >= log2foldchange & padj <= significance), "NLRP3.vs.DEX",   "mmu04621", res.condition.NLRP3.vs.DEX.ann$entrezgene   )




message("+-------------------------------------------------------------------------------")
message("+ UpSetR ")
message("+-------------------------------------------------------------------------------")

library("UpSetR")

listInput <- list(
  'DEX.vs.PLAIN'   = subset(res.condition.DEX.vs.PLAIN.ann,   abs(log2FoldChange) >= log2foldchange & padj <= significance)$ensembl_gene_id, 
  'NLRP3.vs.PLAIN' = subset(res.condition.NLRP3.vs.PLAIN.ann, abs(log2FoldChange) >= log2foldchange & padj <= significance)$ensembl_gene_id, 
  'NLRP3.vs.DEX'   = subset(res.condition.NLRP3.vs.DEX.ann,   abs(log2FoldChange) >= log2foldchange & padj <= significance)$ensembl_gene_id
)


pdf(paste0(Base.dir, "/", Project, "_Fig.UpSetR.pdf"),width=20,height=8, onefile=FALSE)
par(bg=NA)
upset(fromList(listInput),  sets = rev(colnames(fromList(listInput))), nsets = 3, cutoff = 305, # nintersects=NA, 
      keep.order = TRUE, empty.intersections = "on", 
      sets.x.label = paste0("Number of differentially expressed genes\n[l2fc=", log2foldchange, ",adj.p=", significance, "]" ), mainbar.y.label = "Intersections of Differentially Expressed Genes")
dev.off()





message("+-------------------------------------------------------------------------------")
message("+              Deconvolution of the samples using unmix()                       ")
message("+-------------------------------------------------------------------------------")

#Immune_cells_expr_matrix       <- read.csv("srep40508-s1_immuCC_matrix_ensembl.csv")
Immune_cells_expr_matrix            <- read.csv("srep40508-s1.csv")
colnames(Immune_cells_expr_matrix)[1] <- c("external_gene_name")
head(Immune_cells_expr_matrix)

Immune_cells_expr_matrix_anno  <- unique(merge(Immune_cells_expr_matrix, ensEMBL2id, by = "external_gene_name"))
head(Immune_cells_expr_matrix_anno)
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[order(Immune_cells_expr_matrix_anno$entrezgene),]
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[!duplicated(Immune_cells_expr_matrix_anno$external_gene_name),]
rownames(Immune_cells_expr_matrix_anno) <- Immune_cells_expr_matrix_anno$ensembl_gene_id
Immune_cells_expr_matrix_anno  <- Immune_cells_expr_matrix_anno[,-c(1,27:31)]

NormCounts_Immune                 <- as.matrix(normCounts[rownames(normCounts) %in% rownames(Immune_cells_expr_matrix_anno),])
Immune_cells_expr_matrix_anno2    <- (Immune_cells_expr_matrix_anno[rownames(Immune_cells_expr_matrix_anno) %in% rownames(NormCounts_Immune),])

head(NormCounts_Immune)
head(Immune_cells_expr_matrix_anno2)

colnames(NormCounts_Immune)

library(DeconRNASeq)


Decon_results.1_10    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                    signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(1:10)]), 
                                    checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)
Decon_results.1.10_df <- as.data.frame(Decon_results.1_10$out.all)
rownames(Decon_results.1.10_df) <- paste0(sampleTable$condition, "_", sampleTable$sampleName)
#Decon_results.1.10_df$sampleID <- rownames(Decon_results.1.10_df)

Decon_results.11_20    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                    signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(11:20)]), 
                                    checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)
Decon_results.11.20_df <- as.data.frame(Decon_results.11_20$out.all)
rownames(Decon_results.11.20_df) <- paste0(sampleTable$condition, "_", sampleTable$sampleName)
#Decon_results.11.20_df$sampleID <- rownames(Decon_results.11.20_df)


Decon_results.21_25    <- DeconRNASeq(datasets=as.data.frame(NormCounts_Immune), 
                                      signatures=as.data.frame(Immune_cells_expr_matrix_anno2[,c(21:25)]), 
                                      checksig = F, known.prop = FALSE, use.scale = TRUE, fig = TRUE)
Decon_results.21.25_df <- as.data.frame(Decon_results.21_25$out.all)
rownames(Decon_results.21.25_df) <- paste0(sampleTable$condition, "_", sampleTable$sampleName)
#Decon_results.21.25_df$sampleID <- rownames(Decon_results.21.25_df)


Decon_results_df <- merge(Decon_results.1.10_df, Decon_results.11.20_df,by="row.names")
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL

Decon_results_df <- merge(Decon_results_df,  Decon_results.21.25_df,by="row.names",all.x=TRUE)
rownames(Decon_results_df) <- Decon_results_df$Row.names
Decon_results_df$Row.names = NULL


Decon_results_df$condition <- rownames(Decon_results_df)
Decon_results_df$condition <- gsub("_.*", "", Decon_results_df$condition)
head(Decon_results_df)


makeCellTypeDeconvolutionPlot <- function(DECONMATRIX, Cell_Type) {
  t2                <- subset(DECONMATRIX, select=c(Cell_Type, "condition"))
  print(head(t2))
  t2_PLAIN          <- subset(t2, t2$condition == "PLAIN")
  t2_DEX            <- subset(t2, t2$condition == "DEX")
  t2_NLRP3          <- subset(t2, t2$condition == "NLRP3")
  
  colnames(t2_PLAIN)[1] <- "fraction"
  colnames(t2_DEX)[1]   <- "fraction"
  colnames(t2_NLRP3)[1] <- "fraction"
  
  t2.sum_PLAIN <- summarise(group_by(t2_PLAIN, condition), mean=mean(fraction), sd=sd(fraction))
  t2.sum_DEX   <- summarise(group_by(t2_DEX,   condition), mean=mean(fraction), sd=sd(fraction))
  t2.sum_NLRP3 <- summarise(group_by(t2_NLRP3, condition), mean=mean(fraction), sd=sd(fraction))
  
  t2.sum_PLAIN$condition <- "PLAIN"
  t2.sum_DEX$condition   <- "DEX"
  t2.sum_NLRP3$condition <- "NLRP3"
  
 # pdf(paste(Project, "DeconRNASeq_", Cell_Type, "_with_trendline.pdf", sep="_"), width=8, height=6, onefile=FALSE)
#  par(bg=NA)
#  print({
  plot <-   ggplot(t2, aes(x=treatment, y=t2[[1]], fill=condition)) + 
            geom_boxplot(aes(x = condition, y = t2[[1]], group = condition)) + 
      ggtitle(Cell_Type) +
      scale_fill_manual(values = c("purple3",  "#99CCFF", "olivedrab2", "#CCCC33" ,"royalblue2", "#99CCFF", "#339999", "olivedrab2", "#CCCC33" , "slategray","purple3", "royalblue2", "olivedrab2" ,"royalblue2", "#99CCFF", "#339999", "olivedrab2", "#CCCC33" , "slategray") )  + 
      labs(y = "Cell type fraction", x = "condition") +
      geom_line(data=t2.sum_PLAIN, aes(x=condition, y=mean, group = condition), colour="black") + 
      geom_line(data=t2.sum_DEX, aes(x=condition, y=mean, group = condition), colour="black") +
      geom_line(data=t2.sum_NLRP3, aes(x=condition, y=mean, group = condition), colour="black") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
 # })
#  dev.off()
  return(plot)
}


# plt.1 <- makeCellTypeDeconvolutionPlot(Decon_results_df, 'Mast.Cells')
# makeCellTypeDeconvolutionPlot('Neutrophil.Cells')
# makeCellTypeDeconvolutionPlot('Eosinophil.Cells')
# makeCellTypeDeconvolutionPlot('B.Cells.Memory')
# makeCellTypeDeconvolutionPlot('B.Cells.Naive')
# makeCellTypeDeconvolutionPlot('Plasma.Cells')
# makeCellTypeDeconvolutionPlot('T.Cells.CD8.Actived')
# makeCellTypeDeconvolutionPlot('T.Cells.CD8.Naive')
# makeCellTypeDeconvolutionPlot('T.Cells.CD8.Memory')
# makeCellTypeDeconvolutionPlot('M0.Macrophage')
# makeCellTypeDeconvolutionPlot('M1.Macrophage')
# makeCellTypeDeconvolutionPlot('M2.Macrophage')
# makeCellTypeDeconvolutionPlot('Treg.Cells')
# makeCellTypeDeconvolutionPlot('T.Cells.CD4.Memory')
# makeCellTypeDeconvolutionPlot('T.Cells.CD4.Naive')
# makeCellTypeDeconvolutionPlot('T.Cells.CD4.Follicular')
# makeCellTypeDeconvolutionPlot('Th1.Cells')
# makeCellTypeDeconvolutionPlot('Th17.Cells')
# makeCellTypeDeconvolutionPlot('Th2.Cells')
# makeCellTypeDeconvolutionPlot('Monocyte')
# makeCellTypeDeconvolutionPlot('GammaDelta.T.Cells')
# makeCellTypeDeconvolutionPlot('NK.Resting')
# makeCellTypeDeconvolutionPlot('NK.Actived')
# makeCellTypeDeconvolutionPlot('DC.Actived')
# makeCellTypeDeconvolutionPlot('DC.Immature')





message("+-------------------------------------------------------------------------------")
message("+ DeconRNA-Seq : Stacked bar plot                                               ")
message("+-------------------------------------------------------------------------------")

CellTypesList <- c('Plasma.Cells', 'Monocyte', 'Mast.Cells', 'Neutrophil.Cells', 'Eosinophil.Cells',   'NK.Resting', 'NK.Actived', 'M0.Macrophage', 'M1.Macrophage', 'M2.Macrophage','DC.Immature', 'DC.Actived', 'B.Cells.Memory', 'B.Cells.Naive',  'GammaDelta.T.Cells', 'T.Cells.CD8.Actived', 'T.Cells.CD8.Naive', 'T.Cells.CD8.Memory',  'Treg.Cells',  'T.Cells.CD4.Memory', 'T.Cells.CD4.Naive', 'T.Cells.CD4.Follicular', 'Th1.Cells', 'Th17.Cells', 'Th2.Cells'  )
CellTypesList <- gsub("\\.", " ", CellTypesList)

Decon_results_df.means      <- aggregate(Decon_results_df[, c(1:25)], list(Decon_results_df$condition), mean)
Decon_results_df.means.melt <- melt(Decon_results_df.means)

Decon_results_df.vars      <- aggregate(Decon_results_df[, c(1:25)], list(Decon_results_df$condition), var)
Decon_results_df.vars.melt <- melt(Decon_results_df.vars)

Decon_results_df.means.melt$vars     <- Decon_results_df.vars.melt$value

# Calculate a normalised variance
Decon_results_df.means.melt$normvars <- Decon_results_df.means.melt$vars/max(Decon_results_df.means.melt$vars)

# Tidy up the table a bit for plotting
colnames(Decon_results_df.means.melt) <- c("Group", "Immune.Cell", "mean", "variance", "normalisedvariance")
Decon_results_df.means.melt$Immune.Cell <- gsub("\\.", " ", Decon_results_df.means.melt$Immune.Cell)
Decon_results_df.means.melt$Immune.Cell <- factor(Decon_results_df.means.melt$Immune.Cell, levels =  CellTypesList )
head(Decon_results_df.means.melt)



fill <- c("firebrick4", "firebrick1", "tomato2", "orangered", "orange3", "orange1", "gold1", "yellow1", "greenyellow",
          "chartreuse3",  "cyan2", "cyan4","seagreen1","lightskyblue",  "deepskyblue", "dodgerblue", "steelblue2", 
          "royalblue1", "blue3", "navyblue","slateblue2", "slateblue4","darkorchid3", "magenta2", "deeppink2")


Decon.stacked.bar <- ggplot(data = Decon_results_df.means.melt, aes(y = mean, x = Group, fill = Immune.Cell)) + 
                     geom_bar( colour = "black", stat="identity") +  
                     scale_fill_manual(values=fill) +
                     scale_x_discrete(name="", limits=c("PLAIN", "DEX", "NLRP3")) +
                     theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
                     theme(legend.title = element_blank()) + 
                     labs(x="", y="Fraction of cells")
Decon.stacked.bar

pdf(paste(Project, "_DESeq2_ImmuneCellTypeDeconvolution_stacked_barplot.pdf", sep=""), width=10,height=8, onefile=FALSE)
par(bg=NA)
Decon.stacked.bar
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ DeconRNA-Seq : Bubbleplot                                                     ")
message("+-------------------------------------------------------------------------------")

head(Decon_results_df.means.melt)

pdf(paste(Project, "_DESeq2_ImmuneCellTypeDeconvolution_bubbleplot.pdf", sep=""), width=10,height=8, onefile=FALSE)
par(bg=NA)
ggplot(Decon_results_df.means.melt, aes(y = Immune.Cell, x = Group)) +
  geom_point(aes(size = mean, colour = Group, alpha=(1-normalisedvariance))) + 
  xlab("") + ylab("") +
  #scale_alpha_manual(name="Test")  + 
  guides(alpha=guide_legend(title="1-(Normalised Variance)")) +
  guides(size=guide_legend(title="Mean Proportion")) +
  scale_x_discrete(name="", limits=c("PLAIN", "DEX", "NLRP3")) +
  guides(colour=FALSE) +
  scale_color_manual(values = c("PLAIN"="GREY", "DEX"="purple3", "NLRP3"="blue") )  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ Custom Gene Set : Bubbleplots                                                     ")
message("+-------------------------------------------------------------------------------")

makeCustomGeneSetBubbleplots <- function(Base.dir, GENELIST, DDS, ensEMBL2id, Fig.Height, Fig.Width, Fig.ncol, TITLE) {
  
  count            <- 0
  BigTable         <- data.frame(count=double(), condition=character(), gene=character(), genename=character())
  BigTable.summary <- data.frame(condition=character(), mean=double(), gene=character(), genename=character())

  for(gene in GENELIST)  
     {
       message(gene)
       tmp.TC          <- plotCounts(DDS, gene, normalised=TRUE, intgroup = c("condition"), returnData = TRUE)
       tmp.TC$gene     <- gene
       tmp.TC$genename <- subset(ensEMBL2id, ensembl_gene_id == gene)$external_gene_name
  
       tmp.TC.sum          <- ddply(tmp.TC, c("condition"), summarise, mean = median(count))
       tmp.TC.sum$gene     <- gene
       tmp.TC.sum$genename <- subset(ensEMBL2id, ensembl_gene_id == gene)$external_gene_name
  
       BigTable         <- rbind(BigTable, tmp.TC)
       BigTable.summary <- rbind(BigTable.summary, tmp.TC.sum)
       count = count+1 
     }
  
  plt.TC <- ggplot(BigTable, aes(x = condition, y = count, color = condition, group = condition)) + 
            geom_point(data=BigTable.summary, aes(x=as.factor(condition), y=mean, group=condition, colour=condition), size=3, alpha=0.5) + 
            geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
            scale_color_manual(name="Treatment", values = c("PLAIN"="GREY", "DEX"="purple3", "NLRP3"="blue") )  +
            scale_x_discrete(name="", limits=c("PLAIN", "DEX", "NLRP3") )  +
            scale_y_log10() +
            xlab("Treatment") +
            ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
            ggtitle(TITLE) + 
            facet_wrap( ~ genename, ncol=Fig.ncol) +
            guides(size=FALSE) +
            theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                  axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                  strip.background = element_rect(colour="red", fill="whitesmoke"),
                  legend.position="bottom")

  pdf(paste0(Base.dir, "/", Project, "_Fig.TC.", TITLE, ".pdf"),width=Fig.Width,height=Fig.Height)
  par(bg=NA)
  print( plt.TC )
  dev.off()
  
return(plt.TC)
}


plt.1 <- makeCustomGeneSetBubbleplots(Base.dir, Inflammasome_GO_genes_ensembl,          dds, ensEMBL2id, 15, 10, 5, "Inflammasome_GO_genes_ensembl")          # 6
plt.3 <- makeCustomGeneSetBubbleplots(Base.dir, fibrosis_GO_genes_ensembl,              dds, ensEMBL2id, 4,  10, 5, "fibrosis_GO_genes_ensembl")
plt.4 <- makeCustomGeneSetBubbleplots(Base.dir, Fibrosis_myofibroblast,                 dds, ensEMBL2id, 10, 10, 5, "Fibrosis_myofibroblast")
plt.5 <- makeCustomGeneSetBubbleplots(Base.dir, Negative_Regulation_of_Inflammasomes,   dds, ensEMBL2id, 7.5, 10, 5, "Negative_Regulation_of_Inflammasomes")
plt.6 <- makeCustomGeneSetBubbleplots(Base.dir, Signalling_Downstream_of_Inflammasomes, dds, ensEMBL2id, 10, 10, 5, "Signalling_Downstream_of_Inflammasomes") # 4
plt.7 <- makeCustomGeneSetBubbleplots(Base.dir, Main,                                   dds, ensEMBL2id, 7.5, 10, 5, "Main")
plt.8 <- makeCustomGeneSetBubbleplots(Base.dir, Gene.Il10,                              dds, ensEMBL2id,  4, 4,  1, "Gene.Il10")


message("+-------------------------------------------------------------------------------")
message("+ Analysis Complete")
message("+-------------------------------------------------------------------------------")


#------------------------------------------------------------------------------
# FIN
#------------------------------------------------------------------------------

