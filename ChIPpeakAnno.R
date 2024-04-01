#title: Annotation of peaks relative to chromosome regions/features

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ChIPpeakAnno", version = "3.11")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("EnsDb.Hsapiens.v75")
#browseVignettes("ChIPpeakAnno")

## ---- echo=FALSE, results="hide", warning=FALSE--------------------------
suppressPackageStartupMessages({
  library(ChIPpeakAnno)
  library(rtracklayer)
  library(EnsDb.Hsapiens.v75)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(org.Hs.eg.db)
  library(reactome.db)
  library(BSgenome.Hsapiens.UCSC.hg19)
})
knitr::opts_chunk$set(warning=FALSE, message=FALSE)
i <- 0
generateFigureCaption <- function(cap){
  i <<- i+1
  return(paste0("Figure ", i, ". ", cap))
}

## Prepare annotation data with toGRanges
annoData <- toGRanges(EnsDb.Hsapiens.v75, feature="gene")
annoData[1:2]

## Convert peak data to GRanges using toGRages
path <- "~/Fbw7-TF_analysis_Gencode_061220/Hct_Jun_DelWTfcpos_fdr0.01_edgeR_100bp.bed"
peaks <- toGRanges(path, format="BED", header=FALSE) 
peaks[1:2]

## Visualize binding site distribution relative to features 
# Distribution of peaks around TSS
binOverFeature(peaks, annotationData=annoData,
               radius=5000, nbins=20, FUN=length, errFun=0,
               ylab="count", 
               main="Distribution of aggregated peak numbers around TSS")

# Peaks relative to chromosome regions/features

path <- "~/Fbw7-TF_analysis_Gencode_061220/Hct_Jun_DelWTfcpos_fdr0.01_edgeR_100bp.bed"
peaks <- toGRanges(path, format="BED", header=FALSE) 

aCR <- assignChromosomeRegion(peaks, nucleotideLevel=FALSE, 
                            precedence=c("Promoters", "immediateDownstream", 
                                         "fiveUTRs", "threeUTRs", 
                                         "Exons", "Introns"), 
                            TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
barplot(aCR$percentage, las=3)

write.csv(aCR, file = "~/ChIPPeakAnno/Hct_Jun_DelWTfcpos_0.01.csv")

####
peaks <- "~/Intervene/Intervene_HctJun0.01Myc0.05DelWT_fig/sets/0110_Jun Up in Fbw7-Del_Myc Up in Fbw7-Del.bed"
peaks <- toGRanges(path, format="BED", header=FALSE) 

aCR <- assignChromosomeRegion(peaks, nucleotideLevel=FALSE, 
                              precedence=c("Promoters", "immediateDownstream", 
                                           "fiveUTRs", "threeUTRs", 
                                           "Exons", "Introns"), 
                              TxDb=TxDb.Hsapiens.UCSC.hg19.knownGene)
barplot(aCR$percentage, las=3)
