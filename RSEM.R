## Generate a counts file from RSEM data
## The first part is a practice run using data from the tximport vignette
## Scroll for using the tximport package on your own data

##Install Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("tximportData")
BiocManager::install("tximport")

library(tximport)
library(tximportData)
dir <- system.file("extdata", package = "tximportData")
list.files(dir)

samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
samples

files <- file.path(dir, "rsem", samples$run, paste0(samples$run, ".isoforms.results.gz"))
files
names(files) <- paste0("sample", 1:6)
names(files)
txi.rsem <- tximport(files, type = "rsem", txIn = T, txOut = T)
head(txi.rsem$counts)
txi.gene <- summarizeToGene(txi.rsem)

tmpInfo <- data.frame(trt = factor(rep(c("a","b"), each = 3)))
rownames(tmpInfo) <- colnames(txi.rsem$counts)
dds <- DESeqDataSetFromTximport(txi = txi.rsem, colData = tmpInfo, design = ~trt)
vsd <- vst(dds)

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, k, "GENEID", "TXNAME")




###########################################################
## Download the .isoforms files from the NCI PDX Repository
## Then use tximport to process the files and summarize to gene counts
## You can then use DESeqDataSetFromTximport() to process the gene counts for differential expression
## counts() accesses the summarized counts from the dds object which can be used for CMScaller
library(tximport)
library(DESeq2)
##This directory should contain the isoform count files
dir <- list.files("C:/Users/pkuhlers/Documents/rsem")

##This is probably an overly complicated way of telling R the location of all ".isoforms" files
##in the dir
files <- file.path(getwd(), "rsem", dir[grep(".isoforms",list.files("C:/Users/pkuhlers/Documents/rsem"))])
##Make the names whatever you want
names(files) <- paste0("patient", 1:length(files))

##set up tximport object - you want isoform (transcript) counts in and also out (txIn/Out = T)
txi <- tximport(files = files, 
                type = "rsem",
                txIn = T,
                txOut = T)

##This is a table of the HUGO gene names and the corresponding isoforms that make up the gene.
##Its just the first two columns of one of the .isoform files.
##Change the name to one that you have downloaded.
tx2gene <- read.table("rsem/128783_104-T_VJ6DZ2WS2YL1_v2.0.1.4.0_RNASeq.RSEM.isoforms.results",
                       colClasses = c("character", "character", rep("NULL", 6)), sep = "\t", header = T)
txi.gene <- summarizeToGene(txi, tx2gene = tx2gene) ##Adds up all isoform counts to give gene counts

tmpInfo <- data.frame(trt = factor(c("a")))
rownames(tmpInfo) <- colnames(txi.gene$counts)
dds <- DESeqDataSetFromTximport(txi = txi.gene, colData = tmpInfo, design = ~1)

##Classify the samples
##First you need to pull down entrez IDs and match them to the HUGO symbols from the BioMart
library(CMScaller)
cmsRes <- CMScaller(emat = data.frame(counts(dds)),
                    rowNames = "symbol",
                    RNAseq = T,
                    FDR = 0.05)
