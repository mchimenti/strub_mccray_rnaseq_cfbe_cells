## Analysis of McCray/Strubb RNA-seq CFBE cells with delta F508 -/-
## Date: 7.25.2018
## Author: Michael Chimenti
## Organism: hg38 / human
## Aligners: hisat2 / salmon
## Design: Case/control (drug treatment) + known batch effects
## Reps: 4

##########
## Imports
##########

#source("http://bioconductor.org/biocLite.R")

#negative binomial GLM and related
library('DESeq2')
library('calibrate')
library('tximport')
library('readr')
#annotation
library('biomaRt')
library("AnnotationDbi")
library("org.Hs.eg.db")
#Exploratory analysis
library('tidyverse')
library('pcaExplorer')
#pathway
#library(pathview)
#library(gage)
#library(gageData)
#library(ggplot2)

setwd("~/iihg/RNA_seq/mccray/project_strubb_cfbe_july2018/") 

###########
##Function Defs
###########

get_annotation <- function(dds, biomart_dataset, idtype){
  if(is.null(biomart_dataset))
    stop("Select a species to generate the corresponding annotation.
         To obtain a list, type mart = useMart('ensembl'), followed by listDatasets(mart).")
  
  mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="www.ensembl.org",
                  dataset=biomart_dataset)
  anns <- getBM(attributes = c(idtype, "external_gene_name", "description"),
                filters = idtype,
                values = rownames(dds),
                mart = mart)
  
  # keep and match with the ones that are actually there
  anns2 <- anns[match(rownames(dds), anns[, 1]), ]
  rownames(anns2) <- rownames(dds)
  # rename the columns rsp. add row names to be consistent with other function
  colnames(anns2) <- c("gene_id","gene_name","description")
  
  return(anns2)
}

## Volcano Plot function 
volcanoplot <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot", legendpos="topright", labelsig=TRUE, textcx=1, ...) {
  with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=main, ...))
  with(subset(res, padj<sigthresh ), points(log2FoldChange, -log10(pvalue), pch=20, col="red", ...))
  with(subset(res, abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="orange", ...))
  with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), points(log2FoldChange, -log10(pvalue), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, padj<sigthresh & abs(log2FoldChange)>lfcthresh), textxy(log2FoldChange, -log10(pvalue), labs=ext_gene, cex=textcx, offset=0.3, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|LogFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}

#######################################
## tximport > DESeq2 > PCAExplorer
#######################################


samples <- read.table("samples.txt", header=FALSE)
samples$drug <- factor(rep(c("C18","CD1530","DMSO","DORSO","K9499","SALER","VX661","VX809","WITHA","XL147"),each=4))
samples$day <- factor(rep(c("a","b","a","b","b","a","a","a","b","b"),each = 4))  ## day is confounded with drug, so this is not really useful
samples$replicate <- factor(rep(c("rep1","rep2","rep3","rep4"),each = 10))
names(samples) <- c("sample","drug","day")
rownames(samples) <- samples$sample

files <- file.path(getwd(), 'salmon', samples$sample, 'salmon', 'quant.sf')
names(files) <- samples$sample

tx2gene <- read_csv(file.path(getwd(), "salmon", "tx2gene.csv"), col_names = FALSE)

tx2gene$X1 <- tx2gene$X1 %>%
    strsplit(split = '.', fixed = TRUE) %>%
    sapply( "[", 1)  ## obtuse sapply statement needed b/c of annoying way strsplit returns list of lists

txi <- tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ drug)

ddsTxi <- ddsTxi[ rowSums(counts(ddsTxi)) > 5, ]
ddsTxi$drug <- relevel(ddsTxi$drug, ref ='DMSO')

ddsTxi <- DESeq(ddsTxi)

## DE analysis 

res <- results(ddsTxi, contrast = c("drug","WITHA","DMSO"))
res <- na.omit(res)  #drop NA rows
res_sig <- res[res$padj < 0.05 & res$baseMean > 5.0,]
res_ord <- res_sig[order(res_sig$padj),]
res_ord$ext_gene <- anno[row.names(res_ord), "gene_name"]

png("test.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_ord, main = "Volcano Plot:", lfcthresh=1.0, sigthresh=0.1, textcx=.1, xlim=c(-12, 12), ylim = c(0,300))
dev.off()


##---------------launch PCA Explorer on dds object 
anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)

#######
## drop day 2 batches and account for rep 2 batch effect
#######

## drop day 2 batch for PCA reanalysis 
ddsTxi_batch1 <- ddsTxi[,ddsTxi@colData$day == "a"]
ddsTxi_batch1@colData$drug <- droplevels(ddsTxi_batch1@colData$drug)

## add replicate metadata
ddsTxi_batch1@colData$rep <- as.factor(rep(c("one","two","three","four"), times = 5))
ddsTxi_batch1@design <- ~drug + rep
ddsTxi_batch1 <- DESeq(ddsTxi_batch1)

rldTxi_batch1 <- rlog(ddsTxi_batch1, blind = FALSE)
pcaExplorer(dds = ddsTxi_batch1, annotation = anno, rlt = rldTxi_batch1)

res_saler <- results(ddsTxi_batch1, contrast = c("drug", "SALER", "DMSO"))
res_saler <- na.omit(res_saler)
res_saler_sig <- res_saler[res_saler$padj < 0.05 & res_saler$baseMean > 5.0, ]
res_saler_ord <- res_saler_sig[order(res_saler_sig$padj),]
res_saler_ord$ext_gene <- anno[row.names(res_saler_ord), "gene_name"]

res_vx661 <- results(ddsTxi_batch1, contrast = c("drug", "VX661", "DMSO"))
res_vx661 <- na.omit(res_vx661)
res_vx661_sig <- res_vx661[res_vx661$padj < 0.05 & res_vx661$baseMean > 5.0, ]
res_vx661_ord <- res_vx661_sig[order(res_vx661_sig$padj),]
res_vx661_ord$ext_gene <- anno[row.names(res_vx661_ord), "gene_name"]

res_vx809 <- results(ddsTxi_batch1, contrast = c("drug", "VX809", "DMSO"))
res_vx809 <- na.omit(res_vx809)
res_vx809_sig <- res_vx809[res_vx809$padj < 0.05 & res_vx809$baseMean > 5.0, ]
res_vx809_ord <- res_vx809_sig[order(res_vx809_sig$padj),]
res_vx809_ord$ext_gene <- anno[row.names(res_vx809_ord), "gene_name"]

res_c18 <- results(ddsTxi_batch1, contrast = c("drug", "C18", "DMSO"))
res_c18 <- na.omit(res_c18)
res_c18_sig <- res_c18[res_c18$padj < 0.05 & res_c18$baseMean > 5.0, ]
res_c18_ord <- res_c18_sig[order(res_c18_sig$padj),]
res_c18_ord$ext_gene <- anno[row.names(res_c18_ord), "gene_name"]

## Volcano plots 

png("volcano_C18_FDR_0p1.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_c18_ord, main = "Volcano Plot: C18 vs DMSO, FDR p < 0.1", lfcthresh=1.0, sigthresh=0.1, textcx=1, xlim=c(-10, 10), ylim = c(4,20))
dev.off()

png("volcano_vx661_FDR_0p1.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_vx661_ord, main = "Volcano Plot: VX661 vs DMSO, FDR p < 0.1", lfcthresh=1.0, sigthresh=0.1, textcx=0.75, xlim=c(-10, 10), ylim = c(3,40))
dev.off()

png("volcano_vx809_FDR_0p1.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_vx809_ord, main = "Volcano Plot: vx809 vs DMSO FDR < 0.1", lfcthresh=1.0, sigthresh=0.1, textcx=0.75, xlim=c(-10, 10), ylim = c(3,9))
dev.off()

png("volcano_salerimide_FDR_0p1.png", 1200, 1500, pointsize=20, res=100)
volcanoplot(res_saler_ord, main = "Volcano Plot: Salermide vs DMSO, FDR < 0.1", lfcthresh=0.75, sigthresh=0.1, textcx=1, xlim=c(-3, 3), ylim = c(2,15))
dev.off()

## write results 

my_cols <- c("baseMean","log2FoldChange","padj","ext_gene")
write.csv(x = res_c18_ord[,my_cols], file = "DE_genes_C18_DMSO_padj_0p05.csv")
write.csv(x = res_vx809_ord[,my_cols], file = "DE_genes_vx809_padj_0p05.csv")
write.csv(x = res_vx661_ord[,my_cols], file = "DE_genes_vx661_padj_0p05.csv")
write.csv(x = res_saler_ord[,my_cols], file = "DE_genes_salermide_padj_0p05.csv")


