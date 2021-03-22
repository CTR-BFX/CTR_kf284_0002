#!/usr/bin/Rscript

#-------------------------------------------------------------------------------------------------#
# Author  : Dr Russell S. Hamilton
# Email   :	rsh46@cam.ac.uk
# Twitter :	@drrshamilton
# Web     :	http://www.trophoblast.cam.ac.uk/directory/Russell-Hamilton
#
# Project : CTR_kf284_0001
# PI      : Kirstian Franze (kf284@cam.ac.uk)
# Owner   : Damiano Giuseppe Barone (baronedg@gmail.com)
#
# R-Script to perform differential transcript analysis 
#
#-------------------------------------------------------------------------------------------------#
# Copyright Russell S. Hamilton (rsh46@cam.ac.uk) 2017                                            #
# License: 	                                                                                      #
# Attribution-Non Commercial-Share Alike CC BY-NC-SA                                              #
# https://creativecommons.org/licenses/by-nc-sa/                                                  #
#   Attribution:	 You must give appropriate credit, provide a link to the license, and indicate  # 
#                  if changes were made. You may do so in any reasonable manner, but not in any   # 
#                  way that suggests the licensor endorses you or your use.                       #
#   NonCommercial: You may not use the material for commercial purposes.                          #
#   ShareAlike:	   If you remix, transform, or build upon the material, you must distribute your  #
#                  contributions under the same license as the original.                          #
#-------------------------------------------------------------------------------------------------#


# initial install of packages
#source("http://bioconductor.org/biocLite.R")
#biocLite('DESeq2')
#biocLite("limma")

library('DESeq2')
library('RColorBrewer')
library('gplots')
library('ggplot2')
library("ggrepel")
library("RColorBrewer")
library("pheatmap")
library("scales")
library("ggfortify")
library("cowplot")
library("biomaRt")
library("UpSetR")
library("reshape2")
library("VennDiagram")
library("limma")
library("ggalt")
library("dplyr")
library("plyr")
library("BiocParallel")



register(MulticoreParam(3))


elementTextSize  <- 10

l2fc             <- 2
significance     <- 0.01
topN             <- 30


message("+-------------------------------------------------------------------------------")
message("+ Set up some constants e.g. base directories")
message("+-------------------------------------------------------------------------------")

Base.dir <- "~/Documents/CTR-Groups/Kristian_Franze/CTR_kf284_0002_Mouse"
setwd(Base.dir)
HTSeq.dir <- paste(Base.dir,"/HTSEQ-COUNTS", sep="")

Project     <- "CTR_kf284_0002"

message("+-------------------------------------------------------------------------------")
message("+ Set up the sample table")
message("+-------------------------------------------------------------------------------")

sample.table                 <- read.table("sample.table.lanes.txt", header=TRUE)
sample.table$sampleLabelLane <- paste0(sample.table$sampleLabel, "s", sample.table$sampleLane)

sampleTable                  <- data.frame(sampleName=sample.table$sampleLabelLane, 
                                           fileName=sample.table$sampleFiles, 
                                           #lane=sample.table$sampleLabelLane,
                                           label=sample.table$sampleLabel,
                                           condition=sample.table$sampleCondition, 
                                           day=as.factor(sample.table$sampleDay), 
                                           group=as.factor(sample.table$sampleGroup))
str(sampleTable)
nrow(sampleTable)
head(sampleTable)
table( paste0(sampleTable$condition, "_", sampleTable$day))
length(unique( sampleTable$label))


sampleTable$cond_day <- as.factor(paste0(sampleTable$condition, "_d", sampleTable$day) )
head(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Removing samples identified by QC")
message("+-------------------------------------------------------------------------------")
#sampleTable <- sampleTable[-c(3,11),]
nrow(sampleTable)

sampleTable <- subset(sampleTable, label != "G09")

nrow(sampleTable)
str(sampleTable)


message("+-------------------------------------------------------------------------------")
message("+ Retrieve ensEMBL annotations")
message("+-------------------------------------------------------------------------------")
ensembl    =  useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name', 'gene_biotype', 'entrezgene'), mart = ensembl)          
head(ensEMBL2id)

ensEMBL2id_go <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene', 'description', 'go_id', 'name_1006' ), mart = ensembl)    
head(ensEMBL2id_go)

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
message("+ Create ddsHTSeq & dds object")
message("+-------------------------------------------------------------------------------")


ddsHTSeq.condition <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ condition) 
ddsHTSeq.condition <- collapseReplicates(ddsHTSeq.condition, ddsHTSeq.condition$label, renameCols = TRUE)
dds.condition      <- DESeq(ddsHTSeq.condition)


ddsHTSeq.cond_day <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=HTSeq.dir, design= ~ cond_day) 
ddsHTSeq.cond_day <- collapseReplicates(ddsHTSeq.cond_day, ddsHTSeq.cond_day$label, renameCols = TRUE)
dds.cond_day      <- DESeq(ddsHTSeq.cond_day)

message("+-------------------------------------------------------------------------------")
message("+ Inflammasome Plot ")
message("+-------------------------------------------------------------------------------")

Inflammasome.GO.bubble.data <- data.frame("day" = c("d1","d4","d7","d14","d28"), 
                                          "nerve_injury" = c(0.364031,0,0,0,0	), 
                                          "implants" = c(0.319439,0.318472,0.237723,0.333403,0.431433))
Inflammasome.GO.bubble.data.mlt <- melt(Inflammasome.GO.bubble.data)

head(Inflammasome.GO.bubble.data.mlt)

ggplot(Inflammasome.GO.bubble.data.mlt, aes(x = day, y = value, color = variable, group = variable)) + 
  geom_path(aes(x=as.factor(day), y=value, group=variable, colour=variable))  +
  geom_point(aes(x=as.factor(day), y=value, group=variable, colour=variable, size=value), alpha=0.99) + 
  scale_size_area(max_size = 10) +
  scale_x_discrete(limits=c("d1","d4","d7","d14","d28")) +
  scale_color_manual(name="Treatment", values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red")) +
  xlab("Day") +
  ylab(bquote(paste('-log'['10']*'(p-value) GO Enrichment Score'))) +
  guides(size=FALSE) +
  ggtitle("Inflamasome Pathway")

message("+-------------------------------------------------------------------------------")
message("+ DE Analysis ")
message("+-------------------------------------------------------------------------------")

resultsNames(dds.cond_day)
as.data.frame(colData(dds.cond_day))



res.cond_day.nerve_injury_d1.vs.implants_d1   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d1",   "implants_d1"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d4.vs.implants_d4   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d4",   "implants_d4"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d7.vs.implants_d7   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d7",   "implants_d7"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d14.vs.implants_d14 <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d14",  "implants_d14"), type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d28.vs.implants_d28 <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d28",  "implants_d28"), type="normal", parallel=TRUE)

n_nid1_vs_imd1   <- nrow(subset(res.cond_day.nerve_injury_d1.vs.implants_d1,   padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid4_vs_imd4   <- nrow(subset(res.cond_day.nerve_injury_d4.vs.implants_d4,   padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid7_vs_imd7   <- nrow(subset(res.cond_day.nerve_injury_d7.vs.implants_d7,   padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid14_vs_imd14 <- nrow(subset(res.cond_day.nerve_injury_d14.vs.implants_d14, padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid28_vs_imd28 <- nrow(subset(res.cond_day.nerve_injury_d28.vs.implants_d28, padj <= significance & abs(log2FoldChange) >= l2fc) )


res.cond_day.nerve_injury_d1.vs.normal_nerve_d0   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d1",  "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d4.vs.normal_nerve_d0   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d4",  "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d7.vs.normal_nerve_d0   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d7",  "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d14.vs.normal_nerve_d0  <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d14", "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.nerve_injury_d28.vs.normal_nerve_d0  <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "nerve_injury_d28", "normal_nerve_d0"),  type="normal", parallel=TRUE)

n_nid1_vs_nnd0  <- nrow(subset(res.cond_day.nerve_injury_d1.vs.normal_nerve_d0,  padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid4_vs_nnd0  <- nrow(subset(res.cond_day.nerve_injury_d4.vs.normal_nerve_d0,  padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid7_vs_nnd0  <- nrow(subset(res.cond_day.nerve_injury_d7.vs.normal_nerve_d0,  padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid14_vs_nnd0 <- nrow(subset(res.cond_day.nerve_injury_d14.vs.normal_nerve_d0, padj <= significance & abs(log2FoldChange) >= l2fc) )
n_nid28_vs_nnd0 <- nrow(subset(res.cond_day.nerve_injury_d28.vs.normal_nerve_d0, padj <= significance & abs(log2FoldChange) >= l2fc) )


res.cond_day.implants_d1.vs.normal_nerve_d0   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "implants_d1",  "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.implants_d4.vs.normal_nerve_d0   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "implants_d4",  "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.implants_d7.vs.normal_nerve_d0   <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "implants_d7",  "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.implants_d14.vs.normal_nerve_d0  <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "implants_d14", "normal_nerve_d0"),  type="normal", parallel=TRUE)
res.cond_day.implants_d28.vs.normal_nerve_d0  <- lfcShrink(dds=dds.cond_day, contrast=c("cond_day", "implants_d28", "normal_nerve_d0"),  type="normal", parallel=TRUE)

n_imd1_vs_nnd0  <- nrow(subset(res.cond_day.implants_d1.vs.normal_nerve_d0,  padj <= significance & abs(log2FoldChange) >= l2fc) )
n_imd4_vs_nnd0  <- nrow(subset(res.cond_day.implants_d4.vs.normal_nerve_d0,  padj <= significance & abs(log2FoldChange) >= l2fc) )
n_imd7_vs_nnd0  <- nrow(subset(res.cond_day.implants_d7.vs.normal_nerve_d0,  padj <= significance & abs(log2FoldChange) >= l2fc) )
n_imd14_vs_nnd0 <- nrow(subset(res.cond_day.implants_d14.vs.normal_nerve_d0, padj <= significance & abs(log2FoldChange) >= l2fc) )
n_imd28_vs_nnd0 <- nrow(subset(res.cond_day.implants_d28.vs.normal_nerve_d0, padj <= significance & abs(log2FoldChange) >= l2fc) )



DE_Count_Table <- melt( t( data.frame("n_nid1_vs_nnd0"=n_nid1_vs_nnd0, "n_nid4_vs_nnd0"=n_nid4_vs_nnd0, 
                                      "n_nid7_vs_nnd0"=n_nid7_vs_nnd0, "n_nid14_vs_nnd0"=n_nid14_vs_nnd0, 
                                      "n_nid28_vs_nnd0"=n_nid28_vs_nnd0,
                                      "n_imd1_vs_nnd0"=n_imd1_vs_nnd0, "n_imd4_vs_nnd0"=n_imd4_vs_nnd0,
                                      "n_imd7_vs_nnd0"=n_imd7_vs_nnd0, "n_imd14_vs_nnd0"=n_imd14_vs_nnd0,
                                      "n_imd28_vs_nnd0"=n_imd28_vs_nnd0
                                      ) ) )

DE_Count_Table$Var2 <- gsub("n_", "", DE_Count_Table$Var1)
DE_Count_Table$Var2 <- gsub("_vs_.*", "", DE_Count_Table$Var2)
DE_Count_Table$Day  <- gsub(".*d", "", DE_Count_Table$Var2)
DE_Count_Table$Var2 <- gsub("d.*", "", DE_Count_Table$Var2)

DE_Count_Table$Var2 <- gsub("im", "implants", DE_Count_Table$Var2)
DE_Count_Table$Var2 <- gsub("ni", "nerve_injury", DE_Count_Table$Var2)


colnames(DE_Count_Table) <- c("Sample", "Condition", "Counts", "Day")
DE_Count_Table


pdf(paste0(Base.dir, "/", Project, "_Fig.Comparing_Numbers_DEGs_with_time.pdf"),width=5,height=7)
par(bg=NA)
ggplot(data=DE_Count_Table, aes(x=Day, y=Counts, group=Condition, colour=Condition)) +
  geom_point() +
  geom_line() +
  ylab("Number of Differentially Expressed genes \n(Vs normal nerve day 0; adj.pval <= 0.01; lo2FC >=2)") +
#  geom_smooth(aes(group=Condition), method = "glm", se=FALSE) +
  scale_x_discrete(name="Day", limits=c("1","4","7","14","28"))
dev.off()


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



functionPlotDEVolcano <- function(results, title, sig_cut, logfc_cut, topN, col_gender, xrange, yrange) {
  
  results       <- as.data.frame(results)
  results$genes <- results$external_gene_name
  
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log(padj), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    #geom_point(alpha=0.75, size=2, colour=ifelse(results$padj<=sig_cut & results$log2FoldChange >= logfc_cut,"red","grey")) +
    
    geom_point(data=subset(results, abs(results$log2FoldChange) < logfc_cut | results$padj > sig_cut), alpha=0.75, size=1, colour="grey") +
    geom_point(data=subset(results, results$padj<=sig_cut & results$log2FoldChange >= logfc_cut),      alpha=0.75, size=1, colour="red") +
    geom_point(data=subset(results, results$padj<=sig_cut & results$log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=1, colour="blue") +
    
    geom_point(data=subset(results, results$chromosome_name == "Y" & results$padj<=sig_cut & col_gender == 1), alpha=0.99, size=1, colour="blue") +
    geom_point(data=subset(results, results$chromosome_name == "X" & results$padj<=sig_cut & col_gender == 1), alpha=0.99, size=1, colour="pink") +
    
    geom_text_repel( data= subset(results, results$log2FoldChange > 0      & results$padj<=sig_cut)[1:topN,], 
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= tail(subset(results, results$log2FoldChange < 0 & results$padj<=sig_cut),topN),    
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab("log2FC") + ylab(bquote("-log"[10]~"(adj.p.value)")) + 
   scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) + 
  #  scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    #ggtitle(paste0(title, "\n[adj.pval < ", sig_cut, ", abs(l2fc) >", l2fc, ", labelled = top ", topN, "]" )) + 
    ggtitle(paste0(title )) + 
    
    theme(aspect.ratio=1, text = element_text(size=elementTextSize))
  
  return(volc.plt)
}



res.cond_day.nerve_injury_d1.vs.implants_d1.ann   <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d1.vs.implants_d1),   Project, "cond_day.nerve_injury_d1.vs.implants_d1",  significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d4.vs.implants_d4.ann   <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d4.vs.implants_d4),   Project, "cond_day.nerve_injury_d4.vs.implants_d4",  significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d7.vs.implants_d7.ann   <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d7.vs.implants_d7),   Project, "cond_day.nerve_injury_d7.vs.implants_d7",  significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d14.vs.implants_d14.ann <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d14.vs.implants_d14), Project, "cond_day.nerve_injury_d14.vs.implants_d14",significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d28.vs.implants_d28.ann <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d28.vs.implants_d28), Project, "cond_day.nerve_injury_d28.vs.implants_d28",significance,l2fc,ensEMBL2id)


res.cond_day.nerve_injury_d1.vs.implants_d1.ann.vlc   <- functionPlotDEVolcano(res.cond_day.nerve_injury_d1.vs.implants_d1.ann, "nerve_injury_d1 vs implants_d1", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-8,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d4.vs.implants_d4.ann.vlc   <- functionPlotDEVolcano(res.cond_day.nerve_injury_d4.vs.implants_d4.ann, "nerve_injury_d4 vs implants_d4", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-8,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d7.vs.implants_d7.ann.vlc   <- functionPlotDEVolcano(res.cond_day.nerve_injury_d7.vs.implants_d7.ann, "nerve_injury_d7 vs implants_d7", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-8,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d14.vs.implants_d14.ann.vlc <- functionPlotDEVolcano(res.cond_day.nerve_injury_d14.vs.implants_d14.ann, "nerve_injury_d14 vs implants_d14", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-8,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d28.vs.implants_d28.ann.vlc <- functionPlotDEVolcano(res.cond_day.nerve_injury_d28.vs.implants_d28.ann, "nerve_injury_d28 vs implants_d28", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-8,6,2),yrange=c(0,250,50))
pdf(paste0(Base.dir, "/", Project, "_Fig.Volcano_implants-vs-nerve_injury.pdf"),width=5,height=25)
par(bg=NA)
plot_grid(res.cond_day.nerve_injury_d1.vs.implants_d1.ann.vlc, res.cond_day.nerve_injury_d4.vs.implants_d4.ann.vlc, 
          res.cond_day.nerve_injury_d7.vs.implants_d7.ann.vlc, res.cond_day.nerve_injury_d14.vs.implants_d14.ann.vlc, 
          res.cond_day.nerve_injury_d28.vs.implants_d28.ann.vlc, ncol=1)
dev.off()




res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann  <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d1.vs.normal_nerve_d0),  Project, "cond_day.nerve_injury_d1.vs.normal_nerve_d0",  significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d4.vs.normal_nerve_d0.ann  <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d4.vs.normal_nerve_d0),  Project, "cond_day.nerve_injury_d4.vs.normal_nerve_d0",  significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d7.vs.normal_nerve_d0.ann  <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d7.vs.normal_nerve_d0),  Project, "cond_day.nerve_injury_d7.vs.normal_nerve_d0",  significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d14.vs.normal_nerve_d0.ann <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d14.vs.normal_nerve_d0), Project, "cond_day.nerve_injury_d14.vs.normal_nerve_d0", significance,l2fc,ensEMBL2id)
res.cond_day.nerve_injury_d28.vs.normal_nerve_d0.ann <- functionGetSigDEG(as.data.frame(res.cond_day.nerve_injury_d28.vs.normal_nerve_d0), Project, "cond_day.nerve_injury_d28.vs.normal_nerve_d0", significance,l2fc,ensEMBL2id)

res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann.vlc  <- functionPlotDEVolcano(res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann, "nerve_injury_d1 vs normal_nerve_d0", 
                                                                                  significance, l2fc, topN, col_gender=0, xrange=c(-6,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d4.vs.normal_nerve_d0.ann.vlc  <- functionPlotDEVolcano(res.cond_day.nerve_injury_d4.vs.normal_nerve_d0.ann, "nerve_injury_d4 vs normal_nerve_d0", 
                                                                                  significance, l2fc, topN, col_gender=0, xrange=c(-6,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d7.vs.normal_nerve_d0.ann.vlc  <- functionPlotDEVolcano(res.cond_day.nerve_injury_d7.vs.normal_nerve_d0.ann, "nerve_injury_d7 vs normal_nerve_d0", 
                                                                                  significance, l2fc, topN, col_gender=0, xrange=c(-6,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d14.vs.normal_nerve_d0.ann.vlc <- functionPlotDEVolcano(res.cond_day.nerve_injury_d14.vs.normal_nerve_d0.ann, "nerve_injury_d14 vs normal_nerve_d0", 
                                                                                  significance, l2fc, topN, col_gender=0, xrange=c(-6,6,2),yrange=c(0,250,50))
res.cond_day.nerve_injury_d28.vs.normal_nerve_d0.ann.vlc <- functionPlotDEVolcano(res.cond_day.nerve_injury_d28.vs.normal_nerve_d0.ann, "nerve_injury_d28 vs normal_nerve_d0", 
                                                                                  significance, l2fc, topN, col_gender=0, xrange=c(-6,6,2),yrange=c(0,250,50))


pdf(paste0(Base.dir, "/", Project, "_Fig.Volcano_nerve_injury_vs_normal_nerve.pdf"),width=5,height=25)
par(bg=NA)
plot_grid(res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann.vlc, res.cond_day.nerve_injury_d4.vs.normal_nerve_d0.ann.vlc, 
          res.cond_day.nerve_injury_d7.vs.normal_nerve_d0.ann.vlc, res.cond_day.nerve_injury_d14.vs.normal_nerve_d0.ann.vlc, 
          res.cond_day.nerve_injury_d28.vs.normal_nerve_d0.ann.vlc, ncol=1)
dev.off()



res.cond_day.implants_d1.vs.normal_nerve_d0.ann  <- functionGetSigDEG(as.data.frame(res.cond_day.implants_d1.vs.normal_nerve_d0),  Project, "res.cond_day.implants_d1.vs.normal_nerve_d0",  significance,l2fc,ensEMBL2id)
res.cond_day.implants_d4.vs.normal_nerve_d0.ann  <- functionGetSigDEG(as.data.frame(res.cond_day.implants_d4.vs.normal_nerve_d0),  Project, "res.cond_day.implants_d4.vs.normal_nerve_d0",  significance,l2fc,ensEMBL2id)
res.cond_day.implants_d7.vs.normal_nerve_d0.ann  <- functionGetSigDEG(as.data.frame(res.cond_day.implants_d7.vs.normal_nerve_d0),  Project, "res.cond_day.implants_d7.vs.normal_nerve_d0",  significance,l2fc,ensEMBL2id)
res.cond_day.implants_d14.vs.normal_nerve_d0.ann <- functionGetSigDEG(as.data.frame(res.cond_day.implants_d14.vs.normal_nerve_d0), Project, "res.cond_day.implants_d14.vs.normal_nerve_d0", significance,l2fc,ensEMBL2id)
res.cond_day.implants_d28.vs.normal_nerve_d0.ann <- functionGetSigDEG(as.data.frame(res.cond_day.implants_d28.vs.normal_nerve_d0), Project, "res.cond_day.implants_d28.vs.normal_nerve_d0", significance,l2fc,ensEMBL2id)

res.cond_day.implants_d1.vs.normal_nerve_d0.ann.vlc   <- functionPlotDEVolcano(res.cond_day.implants_d1.vs.normal_nerve_d0.ann, "implants_d1 vs normal_nerve_d0", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-6,7.5,2),yrange=c(0,250,50))
res.cond_day.implants_d4.vs.normal_nerve_d0.ann.vlc   <- functionPlotDEVolcano(res.cond_day.implants_d4.vs.normal_nerve_d0.ann, "implants_d4 vs normal_nerve_d0", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-6,7.5,2),yrange=c(0,250,50))
res.cond_day.implants_d7.vs.normal_nerve_d0.ann.vlc   <- functionPlotDEVolcano(res.cond_day.implants_d7.vs.normal_nerve_d0.ann, "implants_d7 vs normal_nerve_d0", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-6,7.5,2),yrange=c(0,250,50))
res.cond_day.implants_d14.vs.normal_nerve_d0.ann.vlc  <- functionPlotDEVolcano(res.cond_day.implants_d14.vs.normal_nerve_d0.ann, "implants_d14 vs normal_nerve_d0", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-6,7.5,2),yrange=c(0,250,50))
res.cond_day.implants_d28.vs.normal_nerve_d0.ann.vlc  <- functionPlotDEVolcano(res.cond_day.implants_d28.vs.normal_nerve_d0.ann, "implants_d28 vs normal_nerve_d0", 
                                                                               significance, l2fc, topN, col_gender=0, xrange=c(-6,7.5,2),yrange=c(0,250,50))

pdf(paste0(Base.dir, "/", Project, "_Fig.Volcano_implants_vs_normal_nerve.pdf"),width=5,height=25)
par(bg=NA)
plot_grid(res.cond_day.implants_d1.vs.normal_nerve_d0.ann.vlc, res.cond_day.implants_d4.vs.normal_nerve_d0.ann.vlc, 
          res.cond_day.implants_d7.vs.normal_nerve_d0.ann.vlc, res.cond_day.implants_d14.vs.normal_nerve_d0.ann.vlc, 
          res.cond_day.implants_d28.vs.normal_nerve_d0.ann.vlc, ncol=1)
dev.off()

message("+-------------------------------------------------------------------------------")
message("+ normCounts")
message("+-------------------------------------------------------------------------------")

normCounts.condition                             <- counts(dds.condition, normalized=TRUE)
normCounts.condition.df                          <- as.data.frame(normCounts.condition)
normCounts.condition.df$ensembl_gene_id          <- rownames(normCounts.condition.df)
normCounts.condition.df.annot                    <- merge(normCounts.condition.df,ensEMBL2id, by="ensembl_gene_id" )
head(normCounts.condition.df.annot)


message("+-------------------------------------------------------------------------------")
message("+ Perform Gender Assignment")
message("+-------------------------------------------------------------------------------")

normCounts.gender           <- normCounts.condition.df.annot
normCounts.gender           <- subset(normCounts.gender, 
                                      (external_gene_name=="Xist" | external_gene_name=="Ddx3y" | 
                                       external_gene_name=="Kdm5d" | external_gene_name=="Uty" | 
                                       external_gene_name=="Zfy1"))

rownames(normCounts.gender) <- normCounts.gender$external_gene_name

normCounts.gender <- normCounts.gender[ , !(names(normCounts.gender) %in% c("ensembl_gene_id","description","entrezgene","gene_biotype"))]
normCounts.gender.mlt <- melt(normCounts.gender)
head(normCounts.gender.mlt)



plt.gender <- ggplot(normCounts.gender.mlt, aes(x=variable, y=log2(value), fill=chromosome_name)) + 
              geom_bar(stat="identity", position="stack") +
              scale_fill_manual(name="Gender Linked Gene Sets", values=(c("X"="pink", "Y"="blue"))) +
              ylab("log2(Normalised Read Counts)") +
              xlab("Samples") +
              theme(text=element_text(size=6,  family="sans"),
                    axis.text.y = element_text(size=6),axis.text.x = element_text(size=6, angle = 90, hjust = 1),
                    legend.position="top")



# plt.gender <- ggplot(normCounts.gender.mlt, aes(x=external_gene_name, y=value, fill=chromosome_name)) + 
#               geom_bar(stat="identity", position=position_dodge()) +
#               coord_flip() +
#               facet_wrap(~ variable, ncol=8) +
#               scale_fill_manual(name="Gender Linked Gene Sets", 
#                     values=(c("X"="pink", "Y"="blue"))) +
#               ylab("Normalised Read Counts") +
#               xlab("Gender Linked Genes") +
#               theme(text=element_text(size=6,  family="sans"),
#                     axis.text.y = element_text(size=6),
#                     axis.text.x = element_text(size=6, angle = 45, hjust = 1),
#                     legend.position="none")


colData(dds.condition)$Gender <- "F"



message("+-------------------------------------------------------------------------------")
message("+ transforms")
message("+-------------------------------------------------------------------------------")

#rld.condition            <- rlogTransformation(dds.condition, blind=T)
vsd.condition            <- varianceStabilizingTransformation(dds.condition, blind=T)
vsd.cond_day             <- varianceStabilizingTransformation(dds.cond_day, blind=T)


message("+-------------------------------------------------------------------------------")
message("+ PCA")
message("+-------------------------------------------------------------------------------")

customPCA <- function(sampleTBL, RLD, TOPNUM, model, ensEMBL2id) {
  
  elementTextSize <- 12
  
  rv     <- rowVars(RLD)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))

  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores <- merge(sampleTBL, pca$x, by="row.names")
  scores$condition_day <- paste0(scores$condition, "_", scores$day)

  scores.summary           <- ddply(scores, c("condition_day"), summarise, meanPC1 = median(PC1), meanPC2 = median(PC2))
  scores.summary$day       <- scores.summary$condition_day
  scores.summary$day       <- gsub(".*_", "",scores.summary$day )
  scores.summary$day       <- as.integer(scores.summary$day)
  scores.summary           <- scores.summary[order(scores.summary$day, decreasing = FALSE),]
  scores.summary$day       <- scores.summary$condition_day
  scores.summary$day       <- gsub(".*_", "",scores.summary$day )
  scores.summary$day       <- as.factor(scores.summary$day)
  scores.summary$condition <- scores.summary$condition_day
  scores.summary$condition <- gsub("_[0-9]+", "",scores.summary$condition )
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, color=day, shape=condition, label=label) ) +
             geom_point(size = 3, alpha=0.75 ) + 
             geom_text_repel(aes(label=label), show.legend = FALSE, size=2, colour="black") +
             scale_shape_manual(name="Treatment", values = c(17, 16, 15)) + 
             xlab(pc1lab) + ylab(pc2lab) + 
             ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
             theme(text = element_text(size=elementTextSize), legend.position="right") +
             guides(condition_day=FALSE) 
    
  scores$day          <- as.factor(scores$day)
  scores.summary$day  <- as.factor(scores.summary$day)
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, shape=condition, group=condition_day) ) +
                geom_encircle(alpha = 0.1, show.legend = FALSE, aes(fill=condition)) +
    
                geom_point(data=scores.summary, aes(x=meanPC1, y=meanPC2, group=condition, colour=condition), shape=5, size=5, alpha=0.75, show.legend=TRUE) + 
                geom_path(data=scores.summary,  aes(x=meanPC1, y=meanPC2, group=condition), show.legend=FALSE)  + 
    
                geom_point(size = 3, alpha=0.75 ) + 
                geom_text_repel(data=scores.summary, aes(x=meanPC1, y=meanPC2, label=condition_day), show.legend = FALSE, size=4, colour="black") +
    
    
                scale_colour_manual(name="Treatment", values=c("normal_nerve"="grey", 
                                                   "nerve_injury"="red",#"nerve_injury_4"="firebrick1","nerve_injury_7"="firebrick3","nerve_injury_14"="red3","nerve_injury_28"="red",
                                                   "implants"="green"#,"implants_4"="lightblue","implants_7"="blue","implants_14"="purple","implants_28"="darkblue")
                                                   )) +
                scale_fill_manual(name="Treatment", values=c("normal_nerve"="grey", "nerve_injury"="red","implants"="green")) +
    
                scale_shape_manual(name="Treatment", values = c(17, 16, 15)) + 
                xlab(pc1lab) + ylab(pc2lab) + #coord_fixed() +
                ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
                theme(text = element_text(size=elementTextSize), legend.position="right") 
  
  loadings                 <- as.data.frame(pca$rotation)
  loadings$ensembl_gene_id <- rownames(loadings)
  loadings                 <- merge(loadings, ensEMBL2id, by="ensembl_gene_id")
  
  pca.1         <- loadings[ order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(pca.1.25$external_gene_name,levels=unique(pca.1.25$external_gene_name)), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[ order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(pca.2.25$external_gene_name,levels=unique(pca.2.25$external_gene_name)), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + theme(axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca, plt.pca.nl, pca.1.25.plot, pca.2.25.plot) )
}

pca.plt.250     <- customPCA(as.data.frame(colData(dds.condition)), assay(vsd.condition), 250,                        "rld.250",  ensEMBL2id)
pca.plt.250[[2]]


pdf(paste0(Base.dir, "/", Project, "_Fig.PCA.Annotatedv2.pdf"),width=10,height=8)
par(bg=NA)
pca.plt.250[[2]]
dev.off()


pca.plt.500     <- customPCA(as.data.frame(colData(dds.condition)), assay(vsd.condition), 500,                        "rld.500",  ensEMBL2id)
pca.plt.1000    <- customPCA(as.data.frame(colData(dds.condition)), assay(vsd.condition), 1000,                       "rld.1000", ensEMBL2id)
pca.plt.all     <- customPCA(as.data.frame(colData(dds.condition)), assay(vsd.condition), nrow(assay(vsd.condition)), "rld.all",  ensEMBL2id)


pdf(paste0(Base.dir, "/", Project, "_Fig.PCA.all_PC1PC2.pdf"),width=15,height=12)
par(bg=NA)
plot_grid(pca.plt.250[[2]], pca.plt.500[[2]], pca.plt.1000[[2]], pca.plt.all[[2]], nrow=2, ncol=2)
dev.off()


pdf(paste0(Base.dir, "/", Project, "_Fig.PCA.500.pdf"),width=10,height=7)
par(bg=NA)
plot_grid(pca.plt.500[[4]], pca.plt.500[[1]], plt.gender, pca.plt.500[[3]], nrow=2, ncol=2, scale = c(0.75, 1, 1, 0.75))
dev.off()







Custom_TimeCourse <- function(dds, time_gene, ensEMBL2id) {

time_gene.ann                <- subset(ensEMBL2id, ensembl_gene_id == time_gene)
timecourse                   <- plotCounts(dds, time_gene, normalised=TRUE, intgroup = c("day","condition"), returnData = TRUE)
timecourse$condition_day     <- paste0(timecourse$condition, "_", timecourse$day)
timecourse                   <- timecourse[order(timecourse$day, decreasing = FALSE),]
timecourse$gene              <- time_gene
timecourse$genename          <- time_gene.ann$external_gene_name

timecourse.summary           <- ddply(timecourse, c("condition_day"), summarise, mean = median(count))
timecourse.summary$day       <- timecourse.summary$condition_day
timecourse.summary$day       <- gsub(".*_", "",timecourse.summary$day )
timecourse.summary$day       <- as.integer(timecourse.summary$day)
timecourse.summary$condition <- timecourse.summary$condition_day
timecourse.summary$condition <- gsub("_[0-9]+", "",timecourse.summary$condition )
timecourse.summary           <- timecourse.summary[order(timecourse.summary$day, decreasing = FALSE),]
timecourse.summary$gene      <- time_gene
timecourse.summary$genename  <- time_gene.ann$external_gene_name

#print(head(timecourse))
#print(head(timecourse.summary))

plt <- ggplot(timecourse, aes(x = day, y = count, color = condition, group = condition)) + 
       geom_path(data=timecourse.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
       geom_point(data=timecourse.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition, size=mean), alpha=0.75) + 
       scale_size_area(max_size = 10) +
       geom_point(size=1, alpha=0.5, shape=1) + 
       scale_color_manual(name="Treatment", values=c("normal_nerve"="black", "nerve_injury"="orange", "implants"="purple"), 
                          limits=c("normal_nerve", "nerve_injury", "implants")) +
       scale_y_log10() +
       xlab("Day") +
       ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
       ggtitle(paste0(time_gene.ann$external_gene_name, "\n[", time_gene, "]")) + 
       guides(size=FALSE) +
       theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
             axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans") )

return(list(plt, timecourse, timecourse.summary))
}



TC.plt.1 <- Custom_TimeCourse(dds.condition, "ENSMUSG00000055717", ensEMBL2id) 
TC.plt.2 <- Custom_TimeCourse(dds.condition, "ENSMUSG00000055717", ensEMBL2id) 

TC.plt.3 <- Custom_TimeCourse(dds.condition, "ENSMUSG00000029379", ensEMBL2id) 
TC.plt.4 <- Custom_TimeCourse(dds.condition, "ENSMUSG00000000204", ensEMBL2id) 
TC.plt.5 <- Custom_TimeCourse(dds.condition, "ENSMUSG00000020897", ensEMBL2id) 
TC.plt.6 <- Custom_TimeCourse(dds.condition, "ENSMUSG00000035385", ensEMBL2id) 


TC.lgnd  <- get_legend(TC.plt.6[[1]]+theme(legend.position="bottom"))

TC.grid  <- plot_grid(TC.plt.3[[1]]+theme(legend.position="none"), TC.plt.4[[1]]+theme(legend.position="none"), 
                      TC.plt.5[[1]]+theme(legend.position="none"), TC.plt.6[[1]]+theme(legend.position="none"), nrow=2, ncol=2)

TC.plt   <- plot_grid( TC.grid, TC.lgnd, ncol=1, rel_heights = c(5, 0.5))


pdf(paste0(Base.dir, "/", Project, "_Fig.TC.grid.pdf"),width=8,height=7)
par(bg=NA)
TC.plt
dev.off()


count            <- 0
BigTable         <- data.frame(count=double(), day=integer(), condition=character(),  
                               condition_day=character(), gene=character(), genename=character())
BigTable.summary <- data.frame(condition_day=character(), mean=double(), day=integer(), 
                               condition=character(), gene=character(), genename=character())

for(gene in Inflammasome_GO_genes_ensembl) 
{
  message(gene)
  tmp.TC <- Custom_TimeCourse(dds.condition, gene, ensEMBL2id) 
  BigTable         <- rbind(BigTable, tmp.TC[[2]])
  BigTable.summary <- rbind(BigTable.summary, tmp.TC[[3]])
  count = count+1 
}


plt.TC.inflamasome <- ggplot(BigTable, aes(x = day, y = count, color = condition, group = condition)) + 
                      geom_path(data=BigTable.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
                      geom_point(data=BigTable.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition), size=3, alpha=0.75) + 
                      geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
                      scale_color_manual(name="Treatment", values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red"), 
                      limits=c("normal_nerve", "nerve_injury", "implants")) +
                      scale_y_log10() +
                      xlab("Day") +
                      ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
                      ggtitle("Inflamasome GO Terms") + 
                      facet_wrap( ~ genename, ncol=5) +
                      guides(size=FALSE) +
                      theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                            axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                            strip.background = element_rect(colour="red", fill="whitesmoke"),
                            legend.position="bottom")
        

pdf(paste0(Base.dir, "/", Project, "_Fig.TC.inflamasome.pdf"),width=10,height=20)
par(bg=NA)
plt.TC.inflamasome
dev.off()



count            <- 0
BigTable         <- data.frame(count=double(), day=integer(), condition=character(),  
                               condition_day=character(), gene=character(), genename=character())
BigTable.summary <- data.frame(condition_day=character(), mean=double(), day=integer(), 
                               condition=character(), gene=character(), genename=character())
Fibrosis_myofibroblast <- c("ENSMUSG00000030968","ENSMUSG00000030970","ENSMUSG00000030983","ENSMUSG00000030986","ENSMUSG00000030987",
                            "ENSMUSG00000031012","ENSMUSG00000031021","ENSMUSG00000080657","ENSMUSG00000080673","ENSMUSG00000035783",	
                            "ENSMUSG00000024563","ENSMUSG00000032402","ENSMUSG00000024515","ENSMUSG00000007613","ENSMUSG00000032440",	
                            "ENSMUSG00000002603")

for(gene in Fibrosis_myofibroblast) 
{
  message(gene)
  tmp.TC <- Custom_TimeCourse(dds.condition, gene, ensEMBL2id) 
  BigTable         <- rbind(BigTable, tmp.TC[[2]])
  BigTable.summary <- rbind(BigTable.summary, tmp.TC[[3]])
  count = count+1 
}


plt.TC.Fibrosis_myofibroblast <- ggplot(BigTable, aes(x = day, y = count, color = condition, group = condition)) + 
                                 geom_path(data=BigTable.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
                                 geom_point(data=BigTable.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition), size=3, alpha=0.75) + 
                                 geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
                                 scale_color_manual(name="Treatment", values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red"), 
                                                    limits=c("normal_nerve", "nerve_injury", "implants")) +
                                 scale_y_log10() +
                                 xlab("Day") +
                                 ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
                                 ggtitle("Fibrosis / myofibroblast genes") + 
                                 facet_wrap( ~ genename, ncol=4) +
                                 guides(size=FALSE) +
                                 theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                                 axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                                 strip.background = element_rect(colour="red", fill="whitesmoke"),
                                 legend.position="bottom")


pdf(paste0(Base.dir, "/", Project, "_Fig.TC.Fibrosis_myofibroblast.pdf"),width=8,height=13.3)
par(bg=NA)
plt.TC.Fibrosis_myofibroblast
dev.off()


#
# Negative_Regulation_of_Inflammasomes
#

count            <- 0
BigTable         <- data.frame(count=double(), day=integer(), condition=character(),  
                               condition_day=character(), gene=character(), genename=character())
BigTable.summary <- data.frame(condition_day=character(), mean=double(), day=integer(), 
                               condition=character(), gene=character(), genename=character())
for(gene in Negative_Regulation_of_Inflammasomes) 
{
  message(gene)
  tmp.TC <- Custom_TimeCourse(dds.condition, gene, ensEMBL2id) 
  BigTable         <- rbind(BigTable, tmp.TC[[2]])
  BigTable.summary <- rbind(BigTable.summary, tmp.TC[[3]])
  count = count+1 
}
plt.TC.NegReg_Inflammasomes <- ggplot(BigTable, aes(x = day, y = count, color = condition, group = condition)) + 
                               geom_path(data=BigTable.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
                               geom_point(data=BigTable.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition), 
                                          size=3, alpha=0.75) + 
                               geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
                               scale_color_manual(name="Treatment", 
                                                  values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red"), 
                                                  limits=c("normal_nerve", "nerve_injury", "implants")) +
                               scale_y_log10() +
                               xlab("Day") +
                               ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
                               ggtitle("Negative Regulation of Inflammasomes") + 
                               facet_wrap( ~ genename, ncol=4) +
                               guides(size=FALSE) +
                               theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                                     axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                                     strip.background = element_rect(colour="red", fill="whitesmoke"),
                                     legend.position="bottom")


pdf(paste0(Base.dir, "/", Project, "_Fig.TC.Negative_Regulation_of_Inflammasomes.pdf"),width=8,height=13.3)
par(bg=NA)
plt.TC.NegReg_Inflammasomes
dev.off()


#
#Signalling_Downstream_of_Inflammasomes
#
count            <- 0
BigTable         <- data.frame(count=double(), day=integer(), condition=character(),  
                               condition_day=character(), gene=character(), genename=character())
BigTable.summary <- data.frame(condition_day=character(), mean=double(), day=integer(), 
                               condition=character(), gene=character(), genename=character())
for(gene in Signalling_Downstream_of_Inflammasomes) 
{
  message(gene)
  tmp.TC <- Custom_TimeCourse(dds.condition, gene, ensEMBL2id) 
  BigTable         <- rbind(BigTable, tmp.TC[[2]])
  BigTable.summary <- rbind(BigTable.summary, tmp.TC[[3]])
  count = count+1 
}
plt.TC.SigDown_Inflammasomes <- ggplot(BigTable, aes(x = day, y = count, color = condition, group = condition)) + 
                                geom_path(data=BigTable.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
                                geom_point(data=BigTable.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition), 
                                           size=3, alpha=0.75) + 
                                geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
                                scale_color_manual(name="Treatment", 
                                                   values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red"), 
                                                   limits=c("normal_nerve", "nerve_injury", "implants")) +
                                scale_y_log10() +
                                xlab("Day") +
                                ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
                                ggtitle("Signalling Downstream of Inflammasomes") + 
                                facet_wrap( ~ genename, ncol=4) +
                                guides(size=FALSE) +
                                theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                                      axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                                      strip.background = element_rect(colour="red", fill="whitesmoke"),
                                      legend.position="bottom")


pdf(paste0(Base.dir, "/", Project, "_Fig.TC.Signalling_Downstream_of_Inflammasomes.pdf"),width=8,height=13.3)
par(bg=NA)
plt.TC.SigDown_Inflammasomes
dev.off()

#
# MAIN
#
count            <- 0
BigTable         <- data.frame(count=double(), day=integer(), condition=character(),  
                               condition_day=character(), gene=character(), genename=character())
BigTable.summary <- data.frame(condition_day=character(), mean=double(), day=integer(), 
                               condition=character(), gene=character(), genename=character())
for(gene in Main) 
{
  message(gene)
  tmp.TC <- Custom_TimeCourse(dds.condition, gene, ensEMBL2id) 
  BigTable         <- rbind(BigTable, tmp.TC[[2]])
  BigTable.summary <- rbind(BigTable.summary, tmp.TC[[3]])
  count = count+1 
}
plt.TC.Main <- ggplot(BigTable, aes(x = day, y = count, color = condition, group = condition)) + 
               geom_path(data=BigTable.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
               geom_point(data=BigTable.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition), 
                          size=3, alpha=0.75) + 
               geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
               scale_color_manual(name="Treatment", 
                                  values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red"), 
                                  limits=c("normal_nerve", "nerve_injury", "implants")) +
               scale_y_log10() +
               xlab("Day") +
               ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
               ggtitle("Main Inflammasomes") + 
               facet_wrap( ~ genename, ncol=4) +
               guides(size=FALSE) +
               theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                     axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                     strip.background = element_rect(colour="red", fill="whitesmoke"),
                    legend.position="bottom")


pdf(paste0(Base.dir, "/", Project, "_Fig.TC.Main.pdf"),width=8,height=13.3)
par(bg=NA)
plt.TC.Main
dev.off()

#
# Gene.Il10
#
count            <- 0
BigTable         <- data.frame(count=double(), day=integer(), condition=character(),  
                               condition_day=character(), gene=character(), genename=character())
BigTable.summary <- data.frame(condition_day=character(), mean=double(), day=integer(), 
                               condition=character(), gene=character(), genename=character())
for(gene in Gene.Il10) 
{
  message(gene)
  tmp.TC <- Custom_TimeCourse(dds.condition, gene, ensEMBL2id) 
  BigTable         <- rbind(BigTable, tmp.TC[[2]])
  BigTable.summary <- rbind(BigTable.summary, tmp.TC[[3]])
  count = count+1 
}
plt.TC.Gene.Il10 <- ggplot(BigTable, aes(x = day, y = count, color = condition, group = condition)) + 
                    geom_path(data=BigTable.summary,  aes(x=as.factor(day), y=mean, group=condition, colour=condition))  +
                    geom_point(data=BigTable.summary, aes(x=as.factor(day), y=mean, group=condition, colour=condition), 
                               size=3, alpha=0.75) + 
                    geom_jitter(size=1, alpha=0.5, shape=1, width=0.1) + 
                    scale_color_manual(name="Treatment", 
                                       values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red"), 
                                       limits=c("normal_nerve", "nerve_injury", "implants")) +
                    scale_y_log10() +
                    xlab("Day") +
                    ylab(bquote(paste('log'['10']*'(Normalised Read Counts)'))) +
                 #   ggtitle("Il10") + 
                    facet_wrap( ~ genename, ncol=1) +
                    guides(size=FALSE) +
                    theme(text=element_text(size=10,family="sans"), plot.title=element_text(size=10,family="sans"),
                          axis.text.x=element_text(size=10,family="sans"), axis.text.y=element_text(size=10,family="sans"),
                          strip.background = element_rect(colour="red", fill="whitesmoke"),
                          legend.position="bottom")


pdf(paste0(Base.dir, "/", Project, "_Fig.TC.Il10.pdf"),width=5,height=4)
par(bg=NA)
plt.TC.Gene.Il10
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ INFLAMASOME IPA BUBBLE PLOT ")
message("+-------------------------------------------------------------------------------")



Inflammasome.GO.bubble.data <- data.frame("day" = c("d1","d4","d7","d14","d28"), 
                                          "nerve_injury" = c(0.364031,0,0,0,0	), 
                                          "implants" = c(0.319439,0.318472,0.237723,0.333403,0.431433))
Inflammasome.GO.bubble.data.mlt <- melt(Inflammasome.GO.bubble.data)

head(Inflammasome.GO.bubble.data.mlt)


pdf(paste0(Base.dir, "/", Project, "_Fig.IPA.inflamasome.bubble pdf"),width=7,height=5)
par(bg=NA)
ggplot(Inflammasome.GO.bubble.data.mlt, aes(x = day, y = value, color = variable, group = variable)) + 
  geom_path(aes(x=as.factor(day), y=value, group=variable, colour=variable))  +
  geom_point(aes(x=as.factor(day), y=value, group=variable, colour=variable, size=2), alpha=0.99) + 
#  scale_size_area(max_size = 10) +
  scale_x_discrete(limits=c("d1","d4","d7","d14","d28")) +
  scale_color_manual(name="Treatment", values=c("normal_nerve"="orange", "nerve_injury"="green", "implants"="red")) +
  xlab("Day") +
  ylab(bquote(paste('-log'['10']*'(p-value) GO Enrighment Score'))) +
  guides(size=FALSE) +
  ggtitle("Inflamasome Pathway")
dev.off()



message("+-------------------------------------------------------------------------------")
message("+ GO Enrichment ")
message("+-------------------------------------------------------------------------------")

library(pathview)
library(gage)
library(gageData)
library(dplyr)
library(cowplot)

source("http://bioconductor.org/biocLite.R")
biocLite('ReactomePA')
#https://bioconductor.org/packages/devel/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
library(ReactomePA)

head(res.cond_day.nerve_injury_d1.vs.implants_d1.ann)

RESULTS <- res.cond_day.nerve_injury_d1.vs.implants_d1.ann

gene_List <- subset(RESULTS, abs(RESULTS$log2FoldChange) > l2fc)
de <- gene_List$entrezgene
head(de)
x <- enrichPathway(gene=de, organism = "mouse", pvalueCutoff=significance, readable=T)
head(as.data.frame(x),20)

png(paste(Project, TITLE , "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff, "barplot.png", sep="_"), width = 1000, height = 600 )
par(bg=NA)
barplot(x, showCategory=8)
dev.off()

png(paste(Project, TITLE , "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff,"dotplot.png", sep="_"), width = 1000, height = 600 )
par(bg=NA)
dotplot(x, showCategory=15)
dev.off()

png(paste(Project, TITLE , "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff,"enrichMap.png", sep="_"), width = 1000, height = 1000 )
par(bg=NA)
cpl <- enrichMap(x, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
dev.off()

png(paste(Project, TITLE , "ReactomePA", "logFoldChanceCutOff", logFoldChanceCutOff,"cnetplot.png", sep="_"), width = 1500, height = 1500 )
par(bg=NA)
cnetplot(x, categorySize="pvalue", foldChange=geneList)
dev.off()

message("+-------------------------------------------------------------------------------")
message("+ KEGG ")
message("+-------------------------------------------------------------------------------")

library(pathview)
library(gage)
library(gageData)
library(dplyr)
library(cowplot)
library("pathview")


data(korg)
head(korg[,1:3], n=20)
data(kegg.sets.mm)
data(go.sets.mm)
data(sigmet.idx.mm)
kegg.sets.mm = kegg.sets.mm[sigmet.idx.mm]
head(kegg.sets.mm, 1)

library(org.Mm.eg.db)

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

FullSet <- res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann$entrezgene

ego <- customKEGG( subset(res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann,  abs(log2FoldChange) >= l2fc & padj <= significance), "nerve_injury_d1.vs.normal_nerve_d0",  "mmu04621", FullSet)
dotplot(ego)


customKEGG( subset(res.cond_day.nerve_injury_d4.vs.normal_nerve_d0.ann,  abs(log2FoldChange) >= l2fc & padj <= significance), "nerve_injury_d4.vs.normal_nerve_d0",  "mmu04621")
customKEGG( subset(res.cond_day.nerve_injury_d7.vs.normal_nerve_d0.ann,  abs(log2FoldChange) >= l2fc & padj <= significance), "nerve_injury_d7.vs.normal_nerve_d0",  "mmu04621")
customKEGG( subset(res.cond_day.nerve_injury_d14.vs.normal_nerve_d0.ann, abs(log2FoldChange) >= l2fc & padj <= significance), "nerve_injury_d14.vs.normal_nerve_d0", "mmu04621")
customKEGG( subset(res.cond_day.nerve_injury_d28.vs.normal_nerve_d0.ann, abs(log2FoldChange) >= l2fc & padj <= significance), "nerve_injury_d28.vs.normal_nerve_d0", "mmu04621")

customKEGG( subset(res.cond_day.implants_d1.vs.normal_nerve_d0.ann,  abs(log2FoldChange) >= l2fc & padj <= significance), "implants_d1.vs.normal_nerve_d0",  "mmu04621")
customKEGG( subset(res.cond_day.implants_d4.vs.normal_nerve_d0.ann,  abs(log2FoldChange) >= l2fc & padj <= significance), "implants_d4.vs.normal_nerve_d0",  "mmu04621")
customKEGG( subset(res.cond_day.implants_d7.vs.normal_nerve_d0.ann,  abs(log2FoldChange) >= l2fc & padj <= significance), "implants_d7.vs.normal_nerve_d0",  "mmu04621")
customKEGG( subset(res.cond_day.implants_d14.vs.normal_nerve_d0.ann, abs(log2FoldChange) >= l2fc & padj <= significance), "implants_d14.vs.normal_nerve_d0", "mmu04621")
customKEGG( subset(res.cond_day.implants_d28.vs.normal_nerve_d0.ann, abs(log2FoldChange) >= l2fc & padj <= significance), "implants_d28.vs.normal_nerve_d0", "mmu04621")





message("+-------------------------------------------------------------------------------")
message("+ UpSetR ")
message("+-------------------------------------------------------------------------------")

library("UpSetR")

listInput <- list(
                  'nerve_injury_d1.vs.normal_nerve_d0'   = subset(res.cond_day.nerve_injury_d1.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'nerve_injury_d4.vs.normal_nerve_d0'   = subset(res.cond_day.nerve_injury_d4.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'nerve_injury_d7.vs.normal_nerve_d0'   = subset(res.cond_day.nerve_injury_d7.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'rnerve_injury_d14.vs.normal_nerve_d0' = subset(res.cond_day.nerve_injury_d14.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'nerve_injury_d28.vs.normal_nerve_d0'  = subset(res.cond_day.nerve_injury_d28.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  
                  'implants_d1.vs.normal_nerve_d0'  = subset(res.cond_day.implants_d1.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'implants_d4.vs.normal_nerve_d0'  = subset(res.cond_day.implants_d4.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'implants_d7.vs.normal_nerve_d0'  = subset(res.cond_day.implants_d7.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'implants_d14.vs.normal_nerve_d0' = subset(res.cond_day.implants_d14.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id, 
                  'implants_d28.vs.normal_nerve_d0' = subset(res.cond_day.implants_d28.vs.normal_nerve_d0.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id,
                  
                  'nerve_injury_d1 vs implants_d1'   = subset(res.cond_day.nerve_injury_d1.vs.implants_d1.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id,
                  'nerve_injury_d4.vs.implants_d4'   = subset(res.cond_day.nerve_injury_d4.vs.implants_d4.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id,
                  'nerve_injury_d7.vs.implants_d7'   = subset(res.cond_day.nerve_injury_d7.vs.implants_d7.ann,   abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id,
                  'nerve_injury_d14.vs.implants_d14' = subset(res.cond_day.nerve_injury_d14.vs.implants_d14.ann, abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id,
                  'nerve_injury_d28.vs.implants_d28' = subset(res.cond_day.nerve_injury_d28.vs.implants_d28.ann, abs(log2FoldChange) >= l2fc & padj <= significance)$ensembl_gene_id
                  )


pdf(paste0(Base.dir, "/", Project, "_Fig.UpSetR.2.pdf"),width=20,height=8, onefile=FALSE)
par(bg=NA)
upset(fromList(listInput),  sets = rev(colnames(fromList(listInput))), nsets = 15, cutoff = 305, # nintersects=NA, 
      keep.order = TRUE, empty.intersections = "on", 
      sets.x.label = paste0("Number of differentially expressed genes\n[l2fc=", l2fc, ",adj.p=", significance, "]" ), mainbar.y.label = "Intersections of Differentially Expressed Genes")
dev.off()


message("+-------------------------------------------------------------------------------")
message("+ END OF SCRIPT ")
message("+-------------------------------------------------------------------------------")
