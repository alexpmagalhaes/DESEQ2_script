################################################################################
### R script for DESeq2                                                       ##
### APMagalhaes                                                               ##
### september 22th, 2015                                                      ##
################################################################################
################################################################################
###                parameters: to be modified by the user                    ###
################################################################################

rm(list=ls())                                        # remove all the objects from the R session

workDir <- "/Users/Alex/Desktop/AthRNAseq_GA_COLD/Results/Results/GA regulated Transcriptome/T21h4VsGAT21h4"      # working directory for the R session

projectName <- "RNAseq_GA_Cold-Ath"                         # name of the project
author <- "APMagalhaes"                                # author of the statistical analysis/report

treated<-"T21h4"                                  # group of interest
untreated<-"GAT21h4"                                  # group to be used as control

################################################################################
###                             running script                               ###
################################################################################
setwd(workDir)
library('DESeq2')

#setup the folder and file structure
directory<-workDir
sampleFiles <- grep("_",list.files(directory),value=TRUE)
sampleFiles

#setup the experimental desing
sampleCondition<-c(treated,treated,treated,untreated,untreated,untreated)
sampleCondition

#load the tables
sampleTable<-data.frame(sampleName=sampleFiles, fileName=sampleFiles, condition=sampleCondition)
sampleTable
#metadata for the experiment
ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design=~condition)
ddsHTSeq
colData(ddsHTSeq)$condition<-factor(colData(ddsHTSeq)$condition, levels=c(untreated,treated))

#DEseq2 analysis
dds<-DESeq(ddsHTSeq)
res<-results(dds)
res<-res[order(res$padj),]
head(res)

#MAPlot
pdf("MAPlot.pdf",width=7,height=7)
plotMA(dds,ylim=c(-2,2),main="DESeq2")
dev.off()

#output DataFrame
mcols(res,use.names=TRUE)
write.csv(as.data.frame(res),file="resultsDESeq2.csv")

#rlog transformation
rld <- rlogTransformation(dds, blind=TRUE)
write.table(as.data.frame(assay(rld),file='DATE-DESeq2-rlog-transformed-counts.txt', sep='\t'))

#heatmap

library("RColorBrewer")
library("gplots")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:30]
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
pdf("rawHeatmap.pdf",width=7,height=7)
heatmap.2(counts(dds,normalized=TRUE)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(20,12))
dev.off()
pdf("rlogHeatmap.pdf",width=7,height=7)
heatmap.2(assay(rld)[select,], col = hmcol,
Rowv = FALSE, Colv = FALSE, scale="none",
dendrogram="none", trace="none", margin=c(20, 12))
dev.off()

#Sample Clustering
pdf("SampleClstr.pdf",width=7,height=7)
distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- colnames(mat) <- with(colData(dds),
paste(condition,sampleFiles , sep=" : "))
hc <- hclust(distsRL)
heatmap.2(mat, Rowv=as.dendrogram(hc),
symm=TRUE, trace="none",
col = rev(hmcol), margin=c(20, 20))
dev.off()

#PCA
pdf("PCAPlot.pdf",width=7,height=7)
print(plotPCA(rld, intgroup=c("condition")))
dev.off()

#Dispertion plot
pdf("DsprtnPlot.pdf",width=7,height=7)
plotDispEsts(dds)
dev.off()

ddsTC <- DESeq(ddsHTSeq, ~ strain + minute + strain:minute)
