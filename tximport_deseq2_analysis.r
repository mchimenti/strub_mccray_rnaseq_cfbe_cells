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

##---------------launch PCA Explorer on dds object 
anno <- get_annotation(ddsTxi, 'hsapiens_gene_ensembl','ensembl_gene_id')
anno <- na.omit(anno)
rldTxi <- rlog(ddsTxi, blind=FALSE)
pcaExplorer(dds=ddsTxi,annotation=anno,rlt=rldTxi)

## look at dispersion estimates 
plotDispEsts(ddsTxi)

ddsTxi_batch1 <- ddsTxi[,ddsTxi@colData$day == "a"]
ddsTxi_batch1@colData$drug <- droplevels(ddsTxi_batch1@colData$drug)
ddsTxi_batch1 <- DESeq(ddsTxi_batch1)

rldTxi_batch1 <- rlog(ddsTxi_batch1, blind = FALSE)
pcaExplorer(dds = ddsTxi_batch1, annotation = anno, rlt = rldTxi_batch1)
