## DeSeq2 practice arena ##

# BiocManager::install("DESeq2")

## Load DESeq2, stringr and ggplot2 libraries

library(DESeq2)
library(plyr)
library(ggplot2)
library(stringr)
library(rebus)
library(readr)
#library(vsn)
#library(gplots)
library(pheatmap)
library(RColorBrewer)
#library(genefilter)
library( org.Hs.eg.db ) 
library(AnnotationDbi) 
#library(EnsDb.Hsapiens.v86) #not loading?
library(GenomicFeatures)
#library(tximport)
library(refGenome)
library(dplyr)
library(ashr)
library(IHW)
library(EnhancedVolcano)
library(biomaRt)
library(ReactomePA)
## Read the output from featurecounts. First set the working directory, then read the output from featurecounts (counts2.txt here) and  save the object.

#setwd("~/example_data/practice_rnaseq_data/featurecounts/")

## Create readcounts and gene transcript/name/ID objects
FT_dir="/home/zacjohnson/OneDrive/Rajagopalan_Lab/RNAseq_new/DEseq_analysis/"
GTF_dir <- "~/OneDrive/Rajagopalan_Lab/RNAseq/raw_data/"

read_counts <- paste(FT_dir,"iHepvPHH_counts.txt",sep="")

cts <- read.delim(read_counts, header=FALSE)
cts_colnames <- unlist(read.table(paste(FT_dir,"sample_names.txt",sep="")))
quick_names <- unlist(read.table(paste(FT_dir,"quick_names.txt",sep="")))
colnames(cts) <- cts_colnames

# annotation attempt
cts$gene_id <- gsub('\\..+$', '', cts$gene_id)

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- cts$gene_id
symbol <- getBM(filters = "ensembl_gene_id",
                attributes = c("ensembl_gene_id","hgnc_symbol"),
                values = genes, 
                mart = mart)
cts <- merge(x = symbol, 
                y = cts, 
                by.x="ensembl_gene_id",
                by.y="gene_id")

cts <- ddply(cts, "hgnc_symbol", numcolwise(sum))
rownames(cts) <- cts$hgnc_symbol
cts <- cts[,-1] #remove ENSGs
cts <- cts[-1,] #remove blanks

test_mat <- data.frame(condition = rep(c("iHep.D0", "iHep.D9", "PHH.D0", "PHH.D9"), each = 4), batch = rep(1,16))[-c(7),]
row.names(test_mat) <- colnames(cts)[-1]
all(rownames(test_mat) %in% colnames(cts))

# Incld outliers DEseq
#dds_cellvtime <- DESeqDataSetFromMatrix(countData = cts, colData = test_mat ,design = ~condition, tidy = TRUE)
#rld <- rlog(dds_cellvtime, blind = FALSE)
#rld_assay <- (assay(rld))
#plotPCA(rld, intgroup="condition")+theme_bw() #looks pretty much the same as the other desingcan

# remove outlier PHH.D0-4
cts_noOutliers <- cts[,!names(cts) %in% c("GRL19940-iHepD0-4-3-7-19_S4")]
tm_noOut <- test_mat[!rownames(test_mat) %in% c("GRL19940-iHepD0-4-3-7-19_S4"),]
all(rownames(tm_noOut) %in% colnames(cts_noOutliers))

## log2-fold (rld) summarization of count data. VST better for large n
dds_cellvtime <- DESeqDataSetFromMatrix(countData = cts_noOutliers, colData = tm_noOut ,design = ~condition, tidy = TRUE)
# remove 0-sum genes 
keep <- rowSums(counts(dds_cellvtime)) > 1
dds_cellvtime <- dds_cellvtime[ keep, ]

# create log-fold permutation of counts
rld <- as.matrix(rlog(dds_cellvtime, blind = FALSE)) #fitType=local for <8 counts, dispersion trend not linear
colnames(rld) <- quick_names
rownames(rld) <- rownames(dds_cellvtime)
PCAs <- prcomp(rld)

#plot PCA of rld to view global differences
DESeq2::plotPCA(rld, intgroup="condition")+theme_bw() #looks pretty much the same as the other design

rld_df <- as.matrix(rld_assay)

# count dispersion summary
dds <- estimateSizeFactors(dds_cellvtime)
df <- bind_rows(
  as_data_frame(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_data_frame(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  


## distribution matrix of sample's Euclidean distance(rlog(expr_geneset))  ##
sampleDistMatrix <- as.matrix(dist(t(assay(rld))))
rownames(sampleDistMatrix) <- rld$condition
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette((brewer.pal(5, "Spectral")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows="correlation",
         clustering_distance_cols="correlation",
         col=colors)

## Clustering of genes within samples (most highly varied genes)
# for pathway enrichment via reactome, top 500 genes
top500Genes <- head(order(-rowVars(rld_df)),500)
mat500 <- rld_df[ top500Genes, ]
mat500 <- mat500 - rowMeans(mat500)
top500comp <- rownames(as.data.frame(prcomp(mat500)$x))
write.table(top500comp,"top500Genes.txt",sep="\t", row.names = F, col.names =F, quote = F)

# Heatmap sorted by variance in genes, 
topVarGenes <- head(order(-rowVars(rld_df)),50)
mat <- rld_df[ topVarGenes, ]
mat <- mat - rowMeans(mat)
pheatmap(mat)

# PCAs of top 50 genes
top50comp <- prcomp(mat)
plot(prcomp(mat))
biplot(prcomp(mat))

PCA_top50 <- as.data.frame(top50comp$x)
topGenes_PCA1 <- PCA_top50[1:10,]
topGenes_PCA2 <- PCA_top50[sort(abs(PCA_top50$PC2),decreasing=T,index.return=T)[[2]],][1:10,]

#Perform LTR to ID genes with consistent profile FROM "Time series experiements" SECTION OF RNA-Seq workflow: gene-level exploratory analysis and differential expression
ddsTC <- DESeq(dds_cellvtime) #LRT to remove factors from "cell_type" and "time_course" which are independent of each other,
res <- results(ddsTC, filterFun=ihw, alpha=0.05, tidy = T)
summary(res)
gene_list <- rownames(rld_df)

plotCounts(ddsTC, gene=which.min(res$padj), intgroup="condition")

# Minimize size effect with ashr in results
resShrink_H9vP9 <- lfcShrink(ddsTC,type = 'ashr', contrast = c("condition","iHep.D9","PHH.D9"))
resShrink_H0vP0 <- lfcShrink(ddsTC,type = 'ashr', contrast = c("condition","iHep.D0","PHH.D0"))
resShrink_H0vH9 <- lfcShrink(ddsTC,type = 'ashr', contrast = c("condition","iHep.D0","iHep.D9"))
resShrink_P0vP9 <- lfcShrink(ddsTC,type = 'ashr', contrast = c("condition","PHH.D0","PHH.D9"))
summary(resShrink_P0vP9)
# plotting for independant filtering, optimal filter = ~20%
plot(metadata(resShrink_P0vP9)$filterNumRej, 
     type="b", ylab="number of rejections",
     xlab="quantiles of filter",
     lines(metadata(resShrink_P0vP9)$lo.fit, col="red"),
     abline(v=metadata(resShrink_P0vP9)$filterTheta))

#plotMA(res.lfc, ylim=c(-10,10), main="LFC test (apeglm)")
## plotting pvalues vs counts
plotMA(resShrink_H9vP9, ylim=c(-10,10), main="iHep D9 v PHH D9")
plotMA(resShrink_H0vP0, ylim=c(-10,10), main="iHep D0 v PHH D0")
plotMA(resShrink_H0vH9, ylim=c(-10,10), main="iHep D0 v iHep D9")
plotMA(resShrink_P0vP9, ylim=c(-10,10), main="PHH D0 v PHH D9")

# could plot pvalues but that's boring
hist(resShrink_P0vP9_test$pvalue[resShrink_P0vP9_test$baseMean > 1], breaks = 0:20/20,
     col = "grey50", border = "white")

# Volcano plots!
EnhancedVolcano(resShrink_P0vP9, lab = rownames(resShrink_P0vP9), x = 'log2FoldChange', y = 'pvalue',
                title='PHH D0 v PHH D9',transcriptPointSize = 2.0,transcriptLabSize = 5.0,
                xlim = c(-5, 8),pCutoff = 0.05, FCcutoff = 1.5)
EnhancedVolcano(resShrink_H0vP0, lab = rownames(resShrink_H0vP0), x = 'log2FoldChange', y = 'pvalue',
                title='iHep D0 v PHH D0',transcriptPointSize = 2.0,transcriptLabSize = 3.0,
                xlim = c(-5, 8),pCutoff = 0.05, FCcutoff = 1.5)
EnhancedVolcano(resShrink_H9vP9, lab = rownames(resShrink_H9vP9), x = 'log2FoldChange', y = 'pvalue',
                title='iHep D9 v PHH D9',transcriptPointSize = 2.0,transcriptLabSize = 3.0,
                xlim = c(-5, 8), ylim = c(-5, 200), pCutoff = 0.05, FCcutoff = 1.5)
EnhancedVolcano(resShrink_H0vH9, lab = rownames(resShrink_H0vH9), x = 'log2FoldChange', y = 'pvalue',
                title='iHep D0 v iHep D9',transcriptPointSize = 2.0,transcriptLabSize = 3.0,
                xlim = c(-5, 8), ylim = c(-5, 200),pCutoff = 0.05, FCcutoff = 1.5)

## export lfc results for reactome (cytoscape) analysis ##
setwd(FT_dir)

# filter for p-value < 0.001 & logfold change > 1.5
lfc_P0vP9 <- na.omit(data.frame(rownames(resShrink_P0vP9),resShrink_P0vP9$log2FoldChange,resShrink_P0vP9$padj))
lfc_P0vP9 <- na.omit(lfc_P0vP9[lfc_P0vP9$resShrink_P0vP9.padj < 0.001 & abs(lfc_P0vP9$resShrink_P0vP9.log2FoldChange) > 1.5,])
lfc_P0vP9 <- lfc_P0vP9[order(lfc_P0vP9$resShrink_P0vP9.log2FoldChange),]
write.table(lfc_P0vP9, file = "P0vP9_allLFC.txt", row.names = F, col.names = F, quote = F)

lfc_P0vH0 <- na.omit(data.frame(rownames(resShrink_H0vP0),resShrink_H0vP0$log2FoldChange,resShrink_H0vP0$padj))
lfc_P0vH0 <- lfc_P0vH0[lfc_P0vH0$resShrink_H0vP0.padj < 0.001 & abs(lfc_P0vH0$resShrink_H0vP0.log2FoldChange) > 1.5,]
lfc_P0vH0 <- lfc_P0vH0[order(lfc_P0vH0$resShrink_H0vP0.log2FoldChange),]
write.table(lfc_P0vH0, file = "P0vH0_allLFC.txt", row.names = F, col.names = F, quote = F)

lfc_H9vP9 <- data.frame(rownames(resShrink_H9vP9),resShrink_H9vP9$log2FoldChange,resShrink_H9vP9$padj)
lfc_H9vP9 <- na.omit(lfc_H9vP9[lfc_H9vP9$resShrink_H9vP9.padj < 0.001 & abs(lfc_H9vP9$resShrink_H9vP9.log2FoldChange) > 1.5,])
lfc_H9vP9 <- lfc_H9vP9[order(lfc_H9vP9$resShrink_H9vP9.log2FoldChange),]
write.table(lfc_H9vP9, file = "H9vP9_allLFC.txt", row.names = F, col.names = F, quote = F)

lfc_H0vH9 <- data.frame(rownames(resShrink_H0vH9),resShrink_H0vH9$log2FoldChange,resShrink_H0vP0$padj)
lfc_H0vH9 <- na.omit(lfc_H0vH9[lfc_H0vH9$resShrink_H0vP0.padj < 0.001 & abs(lfc_H0vH9$resShrink_H0vH9.log2FoldChange) > 1.5,])
lfc_H0vH9 <- lfc_H0vH9[order(lfc_H0vH9$resShrink_H0vH9.log2FoldChange),]
write.table(lfc_H0vH9, file = "H0vH9_allLFC.txt", row.names = F, col.names = F, quote = F)

# Up & Down regulated genes
lfc_up_P0vP9 <- lfc_P0vP9[lfc_P0vP9$resShrink_P0vP9.log2FoldChange > 0,]
write.table(lfc_up_P0vP9, file = "P0vP9_up_LFC.txt", row.names = F, col.names = F, quote = F)
lfc_down_P0vP9 <- lfc_P0vP9[lfc_P0vP9$resShrink_P0vP9.log2FoldChange < 0,]
write.table(lfc_down_P0vP9, file = "P0vP9_down_LFC.txt", row.names = F, col.names = F, quote = F)

lfc_up_P0vH0 <- lfc_P0vH0[lfc_P0vH0$resShrink_H0vP0.log2FoldChange > 0,]
write.table(lfc_up_P0vH0, file = "H0vP0_up_LFC.txt", row.names = F, col.names = F, quote = F)
lfc_down_P0vH0 <- lfc_P0vH0[lfc_P0vH0$resShrink_H0vP0.log2FoldChange < 0,]
write.table(lfc_down_P0vH0, file = "H0vP0_down_LFC.txt", row.names = F, col.names = F, quote = F)

lfc_up_H9vP9 <- lfc_H9vP9[lfc_H9vP9$resShrink_H9vP9.log2FoldChange > 0,]
write.table(lfc_up_H0vH9, file = "H9vP9_up_LFC.txt", row.names = F, col.names = F, quote = F)
lfc_down_H9vP9 <- lfc_H9vP9[lfc_H9vP9$resShrink_H9vP9.log2FoldChange < 0,]
write.table(lfc_down_H9vP9, file = "H9vP9_down_LFC.txt", row.names = F, col.names = F, quote = F)

lfc_up_H0vH9 <- lfc_H0vH9[lfc_H0vH9$resShrink_H0vH9.log2FoldChange > 0,]
write.table(lfc_up_H0vH9, file = "H0vH9_up_LFC.txt", row.names = F, col.names = F, quote = F)
lfc_down_H0vH9 <- lfc_H0vH9[lfc_H0vH9$resShrink_H0vH9.log2FoldChange < 0,]
write.table(lfc_down_H0vH9, file = "H0vH9_down_LFC.txt", row.names = F, col.names = F, quote = F)

