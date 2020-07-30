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


#4.1
nrow(dds)
eep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
# at least 3 samples with a count of 10 or higher
keep <- rowSums(counts(dds) >= 10) >= 3

#5.1
dds <- DESeq2(dds)

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

#6.1
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("dex"))
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("dex","cell"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = dex, y = count, color = cell)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)
ggplot(geneCounts, aes(x = dex, y = count, color = cell, group = cell)) +
  scale_y_log10() + geom_point(size = 3) + geom_line()

