#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#  BiocManager::install("DESeq2", version = "3.16")

library("DESeq2")
dir <- "C:/Users/anden/Desktop/X4/P4/GenAn/Lab/htseq/"

#### Differential expression between species ####

sampleFiles <- grep("Aril",list.files(dir),value=TRUE)
condition <- c('Monthong','Monthong','Monthong','MusangK','MusangK')

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir,
                                       design= ~ condition)
# See data object
ddsHTSeq

# Filtering
dds <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]
summary(res)

plotMA(res, main="DESeq2", ylim=c(-2,2))

#install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:10]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[c("condition")])
labels <- c("E3 ubiquitin-protein ligase", "sialylation","V-type proton ATPase subunit", "gag-polypeptide of LTR copia-type", "UNKNOWN", "methionine", "glutamate synthase", "UNKNOWN", "UNKNOWN", "xyloglucan galactosyltransferase")
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, labels_row=labels)

rld <- rlog(dds)
plotPCA(rld, intgroup=c("condition"))
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))

#### Differential expression between plant organs ####

sampleFiles <- grep("MusangK",list.files(dir),value=TRUE)
condition <- c('Aril','Aril','Feaf','Root','Stem')

sampleTable <- data.frame(sampleName = sampleFiles,
                          fileName = sampleFiles,
                          condition = condition)

ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                       directory = dir,
                                       design= ~ condition)
# See data object
ddsHTSeq

# Filtering
dds <- ddsHTSeq[ rowSums(counts(ddsHTSeq)) > 1, ]

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$padj),]
summary(res)

plotMA(res, main="DESeq2", ylim=c(-2,2))

#install.packages("pheatmap")
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:10]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[c("condition")])
labels <- c("E3 ubiquitin-protein ligase", "UNKNOWN", "V-type proton ATPase subunit", "UNKNOWN", "methionine", "UNKNOWN", "UNKNOWN", "glutamate synthase", "xyloglucan galactosyltransferase", "UNKNOWN")
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, labels_row=labels)


rld <- rlog(dds)
plotPCA(rld, intgroup=c("condition"))
data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=condition, shape=type)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
#### End ####

