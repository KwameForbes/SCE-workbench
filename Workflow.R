
#2.3
library("airway")
dir <-system.file("extdata", package = "airway", mustWork = TRUE )
list.files(dir)
list.files(file.path(dir,"quants"))
csvfile <- file.path(dir, "sample_table.csv")
coldata <- read.csv(csvfile, row.names = 1, stringsAsFactors = FALSE)
coldata
coldata <- coldata[1:2,]
coldata$names <- coldata$Run
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf.gz")
file.exists(coldata$files)
library("tximeta")
se <- tximeta(coldata)
dim(se)
head(rownames(se))

#2.5
data(gse)
gse

#3
gse$donor
gse$condition
gse$cell <- gse$donor
gse$dex <- gse$condition
levels(gse$dex)
levels(gse$dex) <- c("untrt", "trt")
library("magrittr")
gse$dex %<>% relevel("untrt")
gse$dex
gse$dex <- relevel(gse$dex, "untrt")

#3.1
library("DESeq2")
dds <- DESeqDataSet(gse, design = ~ cell + dex)

#4.1
nrow(dds)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 3

#5.1
dds <- DESeq(dds)

#5.2
res <- results(dds)
res <- results(dds, contrast=c("dex","trt","untrt"))
mcols(res, use.names = TRUE)
summary(res)
res.05 <- results(dds, alpha = 0.05)
table(res.05$padj < 0.05)
resLFC1 <- results(dds, lfcThreshold=1)
table(resLFC1$padj < 0.1)

#5.3
results(dds, contrast = c("cell", "N061011", "N61311"))

#5.4
sum(res$pvalue < 0.05, na.rm=TRUE)
sum(!is.na(res$pvalue))
sum(res$padj < 0.1, na.rm=TRUE)
resSig <- subset(res, padj < 0.1)
head(resSig[ order(resSig$log2FoldChange), ])
head(resSig[ order(resSig$log2FoldChange, decreasing = TRUE), ])
head(rownames(res))
rownames(res) <- sub("\\..*","",rownames(res))




#6.1
topGene <- rownames(res)[which.min(res$padj)]
rownames(dds) <- rownames(res)
plotCounts(dds, gene = my.gene, intgroup=c("dex"))
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

#Take e.g. topGene and then see if you can match up the gene with a single cell dataset.
#so because airway is human, we will need to find a human single cell dataset. we can try the PBMC dataset

BiocManager::install("TENxPBMCData")
library(TENxPBMCData)
tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
tenx_pbmc4k
args(TENxPBMCData)
counts(tenx_pbmc4k)

options(DelayedArray.auto.block.size = 1e9)
lib.sizes <- colSums(counts(tenx_pbmc4k))
n.exprs <- colSums(counts(tenx_pbmc4k) != 0L)
ave.exprs <- rowMeans(counts(tenx_pbmc4k))

destination <- tempfile()
saveRDS(tenx_pbmc4k, file = destination)

sessionInfo()
sce <- tenx_pbmc4k
sce <- logNormCounts(sce)

sce <- runPCA(sce)
sce <- runTSNE(sce)
sce <- runUMAP(sce)

library(scran)
g <- buildSNNGraph(sce, k=10, use.dimred = 'PCA')
clust <- igraph::cluster_walktrap(g)$membership
table(clust)

library(scater)
colLabels(sce) <- pred$labels
#plotReducedDim(sce, "TSNE", colour_by="label")


#library(scater)
min_adj_pval <- which.min(res$padj)
min_adj_pval
rownames(res)[min_adj_pval]
#library(SingleCellExperiment)

#my.gene <- rownames(res)[min_adj_pval]
o <- order(res$padj)
my.gene <- rownames(res)[o[3]]
my.gene %in% rownames(sce)
sum(counts(sce[my.gene,]))


logcts <- logcounts(sce)[my.gene,]
plotColData(sce, y=I(logcts), x="label")

plotCounts(dds, gene = my.gene, intgroup=c("dex"))


#Annotating
BiocManager::install("SingleR")
BiocManager::install("org.Hs.eg.db")
BiocManager::install('pheatmap')

sce2 <- sce
library(org.Hs.eg.db)
rownames(sce2) <- mapIds(org.Hs.eg.db, rownames(sce2), "SYMBOL", "ENSEMBL")

library(SingleR)

ref <- BlueprintEncodeData()
pred <- SingleR(test=sce2, ref=ref, labels=ref$label.main)
table(pred$labels)

table(pred$labels, colLabels(sce))

plotScoreHeatmap(pred)

integrateWithSingleCell<- function(res, dds) {
  # figure out organism from dds
  s <- rownames(dds)
  p <- startsWith(s, "ENSG")
  r <- startsWith(s, "ENSMUSG")
  if ( p[1] == TRUE) {
    print("Your dataset appears to be Human.")
    
    
  }else if (r[1] == TRUE) {
    print("Your dataset appears to be mouse.")
    
  }else {
    print("We only support human and mouse datasets.")
  }  
  tab1 <- data.frame(name=c("pbmc4k","pbmc8k"),
                     pub=c("Hansen 2020","Hansen 2020"),
                     nCells=c(4340,8381),
                     description=c("PBMCs","PBMCs"))
  tab2 <- data.frame(name=c("testing","testing"),
                     pub=c("Hansen 2020","Hansen 2020"),
                     nCells=c(0,0),
                     description=c("PBMCs","PBMCs"))
  pbmc4k <- function(){
    #BiocManager::install("TENxPBMCData")
    #library(TENxPBMCData)
    tenx_pbmc4k <- TENxPBMCData(dataset = "pbmc4k")
    sce <- tenx_pbmc4k
    args(TENxPBMCData)
    counts(tenx_pbmc4k)
    return(sce)
  }
  if (p[1] == TRUE) {
    print("Choose a human single-cell to integrate with your dataset.")
    print(tab1)
    ans <- menu(tab1$name)
    sce <- do.call(tab1$name[ans], list())
  }else if (r[1]== TRUE) {
    print("Choose a mouse single-cell to integrate with your dataset.")
    print(tab2)
    ans <- menu(tab2$name)
    sce <- do.call(tab2$name[ans], list())
  }
  sce
  #sce <- do.call(tab$name[ans], list())
  #print(tab)
  #ans <- menu(tab$name)
  
  #return(list(res=res, dds=dds, ans=ans))
}
# provide relevant single cell dataset to user
# do all your hard work
#return(list(res=res, dds=dds, ans=ans))
plotter <- function(dat) {
  stopifnot(all(names(dat) == c("res", "dds", "ans")))
  plot(dat$res, dat$dds)
}
# example code:
dat <- foo(res, dds)
plotter(dat)
