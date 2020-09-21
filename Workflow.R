#WORKFLOW
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
#head(rownames(se))

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
  
  org <- if (startsWith(rownames(dds[1]), "ENSG")) "human" else "mouse" 
  
  print(paste("Your dataset appears to be",org))
  
  tab <- data.frame(
    pkg=c("TENxPBMCData","TENxPBMCData") ### TODO
    func=c("BaronPancreasData","BaronPancreasData",,"test"), ### TODO
    data=c("human","mouse","pbmc4k","pbmc8k","human"),
    orgm=c("human","mouse","human","human","human"),
    pub=c("Hansen 2020","Hansen 2020","Hansen 2020","Hansen 2020","test 2020"),
    nCells=c(4340,8381,4340,8381,0000),
    desc=c("PBMCs","PBMCs","PBMCs","PBMCs","Red")
  )
  
  print(paste("Choose a",org,"single-cell dataset to integrate with."))
  tab <- tab[tab$orgm == org,]
  tab2 <- tab[,c("pkg","func","data", "pub","nCells","desc")]
  rownames(tab2) <- seq_len(nrow(tab2))
  print(tab2)
  ans <-menu(paste(tab2$func, tab2$data, sep="-"))
  pkg <- tab$pkg[ans]
  if (!requireNamespace(package=pkg, quietly=TRUE)) {
    message(paste0("Package: '",pkg, "' not installed."))
    ask <- askYesNo("Would you like to install package?")
    if (ask) {
      BiocManager::install(pkg)
    } else {
      stop("Package would need to be installed.")
    }
    if (!requireNamespace(package=pkg, quietly=TRUE)) {
      message("Package installed successfully.")
    } else {
      stop("Package needs to be installed.")
    }
  }

  # load package
  require(pkg, character.only = TRUE)
 

  # if the dataset is in the scRNAseq package...
  if (pkg == "scRNAseq") {
    # if only one dataset within the function...
    if (is.na(tab$data[ans])) {
      print("HI")
      sce <- do.call(tab$func[ans], list(ensembl=TRUE))
    } else {
      print("HI2")
      sce <- do.call(tab$func[ans], list(which=tab$data[ans], ensembl=TRUE))
    }
  } else {
    # if only one dataset within the function...
    if (is.na(tab$data[ans])) {
      print("HI3")
      sce <- do.call(tab$func[ans], list())
    } else {
      print("HI4")
      sce <- do.call(tab$func[ans], list(dataset=tab$data[ans]))
    }
  }
  return(list(res=res, dds=dds, sce=sce))
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

