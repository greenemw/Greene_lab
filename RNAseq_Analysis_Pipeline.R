##This is an analysis pipeline for RNA-seq data using the DESeq2 package.
##It will take a RNA-seq count matrix and evaluate differential expression, principal components
##and prepare data for export and further analysis.

##Install Bioconductor
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("TCGAbiolinks")

library(DESeq2)
library(TCGAbiolinks)
library(tidyverse)
library(ggplot2)


##Choose the files that contains the RNAseq counts and the associated phenotype table (info table)
PDXCOUNTSFILE <- file.choose()
PDXINFOFILE <- file.choose()

##Read in the RNA-seq Counts and infotable
PDXCOUNTS <- read.table(PDXCOUNTSFILE,header = T,row.names = 1)
PDXINFOTABLE <- read.csv(PDXINFOFILE)

##Make sure that the names in the info table match up with the column names of the counts table
PDXCOUNTS <- PDXCOUNTS[, as.character(PDXINFOTABLE$Description)]
all(colnames(PDXCOUNTS) == PDXINFOTABLE$Description) ##Sanity check


##Its easier to do contrasts in the PDX line and tissue type is combined, so add that to the info table
PDXINFOTABLE$Line_Type <- paste(PDXINFOTABLE$PDX_Line, PDXINFOTABLE$Tumor_type, sep = "_")
##Make sure that column is a factor since its going into DESeq
PDXINFOTABLE$Line_Type <- factor(PDXINFOTABLE$Line_Type)

##Now set up the dds object
dds <- DESeqDataSetFromMatrix(countData = PDXCOUNTS,
                              colData = PDXINFOTABLE,
                              design = ~Line_Type)
##Run DESeq
dds <- DESeq(dds)

##Use results to make different contrasts between the groups in Line_Type
res1003 <- results(dds, contrast = c("Line_Type", "1003_3D", "1003_SCID"))



##Annotating Genes in the results matrices
##Remember to change dataset and attributes to the species you are working with
##Install from Bioconductor if necessary
#BiocManager::install("biomaRt")
library(biomartr)
library(biomaRt)

##This will list the available biomaRt datasets
# mart <- useEnsembl("ENSEMBL_MART_ENSEMBL")
# listDatasets(mart)

##This will list what attributes you can pull from the respective mart
# listAttributes(mart)

##Human biomart
mart <- useEnsembl(biomart = "ensembl", 
                   dataset = "hsapiens_gene_ensembl")
bm <- getBM(attribute = c("ensembl_gene_id", "hgnc_symbol", "description"), 
            mart = mart)

##Mouse biomart
#mart <- useEnsembl(biomart = "ensembl",
#                   dataset = "mmusculus_gene_ensembl" )
#bm <- getBM(attribute = c("ensembl_gene_id", "mgi_symbol", "mgi_description"), 
#            mart = mart)


##This will merge the DESeq results matrix with annotated gene data. Do this for all result tables.
##This will also create a csv of the results table to be used for other analysis
res1003$ensembl_gene_id <- rownames(res1003)
geneRes1003 <- merge(data.frame(res1003), bm)
write.csv(geneRes1003, "File_Name.csv")


##Exploratory Analyses - use the variance stabilized data
vsd <- vst(dds, blind = F)

plotPCA(vsd, intgroup = "Line_Type")



#########################################################################
##Clasification of sample subtypes with CMScaller
library(devtools)
install_github("LotheLab/CMScaller")
library(CMScaller)

##Use CMScaller to subtype samples
callerRes <- CMScaller(emat = PDXCOUNTS, 
                       rowNames = "ensg",
                       RNAseq = T, 
                       FDR = 0.05)

##Perform gene set analysis with the subtyped samples, first convert ensemble to entrez.
entrezCount <- replaceGeneId(emat = PDXCOUNTS,
                             id.in = "ensg",
                             id.out = "entrez")

gsa <- CMSgsa(emat = entrezCount, 
              class = callerRes$prediction,
              RNAseq = T)

#############################################################################
##Compare the PDX data with the TCGA data - DEGs, PCA

##Queries HT-Seq counts data from COAD project
query.RNAseq.COAD <- GDCquery(project = "TCGA-COAD", 
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "HTSeq - Counts")
##Downloads the query above
GDCdownload(query.RNAseq.COAD)
##Puts all of the data together into a summarized experiment object, adds in the clinical data as well.
COADRNAseqData <- GDCprepare(query = query.RNAseq.COAD,
                             summarizedExperiment = TRUE)

##Put the TCGA counts data into a data frame that we can later merge with the PDX data
TCGACOUNTS <- data.frame(assay(COADRNAseqData))
colnames(TCGACOUNTS) <- colnames(assay(COADRNAseqData))

##Merge TCGA with PDX by adding an ensembl column to both
pdxCounts <- PDXCOUNTS
pdxCounts$ensembl <- rownames(PDXCOUNTS)  
tcgaCounts <- TCGACOUNTS
tcgaCounts$ensembl <- rownames(TCGACOUNTS)

pdx_tcga_counts <- merge(pdxCounts,tcgaCounts, by = "ensembl")
pdx_tcga_counts <- column_to_rownames(pdx_tcga_counts, var = "ensembl")
colnames(pdx_tcga_counts) <- gsub(pattern = "-", replacement = "_", x = colnames(pdx_tcga_counts))

##Now construct an info table with relevant information from TCGA and data in the PDX info table
pdx_tcga_names <- c(as.character(PDXINFOTABLE$Description), COADRNAseqData$barcode)
pdx_tcga_names <- gsub(pattern = "-",
                       replacement = "_",
                       pdx_tcga_names)

pdx_tcga_type <- c(as.character(PDXINFOTABLE$Line_Type), COADRNAseqData$definition)
pdx_tcga_type <- gsub(pattern = " ", 
                      replacement = "_",
                      x = pdx_tcga_type) 

PDX_TCGA_INFOTABLE <- data.frame(name = pdx_tcga_names, type = pdx_tcga_type) 
##Need to remove 'metastatic' and 'recurrent' from the info table - only 1 sample each
PDX_TCGA_INFOTABLE <- subset(x = PDX_TCGA_INFOTABLE,subset = type !="Recurrent_Solid_Tumor" & 
                               type != "Metastatic")

PDX_TCGA_COUNTS <- pdx_tcga_counts[,as.character(PDX_TCGA_INFOTABLE$name)]
all(colnames(PDX_TCGA_COUNTS) == PDX_TCGA_INFOTABLE$name) ##sanity check
#which(colnames(PDX_TCGA_COUNTS) != PDX_TCGA_INFOTABLE$name) ##sanity check

pdx_tcga_dds <- DESeqDataSetFromMatrix(countData = PDX_TCGA_COUNTS,
                                       colData = PDX_TCGA_INFOTABLE,
                                       design = ~type)
pdx_tcga_dds <- DESeq(pdx_tcga_dds)



res_pdx3D1003_tcgaNormal <- results(pdx_tcga_dds, 
                                    contrast = c("type", "1003_3D", "Solid_Tissue_Normal")) 


##Annotate the results table with gene information
res_pdx3D1003_tcgaNormal$ensembl_gene_id <- rownames(res_pdx3D1003_tcgaNormal)
geneRes_pdx3D1003_tcgaNormal <- merge(data.frame(res_pdx3D1003_tcgaNormal), bm, by = "ensembl_gene_id")




##Variance Stabilize for downstream analysis
pdx_tcga_vsd <- vst(pdx_tcga_dds)

plotPCA(pdx_tcga_vsd, intgroup = c("type"))
##Use the z-score function on VST data for GSEA analysis
zScore <- function(gene_counts){
  rMean <- rowMeans2(gene_counts)
  rStDev <- rowSds(gene_counts)
  zScore <- (gene_counts - rMean) / rStDev
  return(zScore)
}

pdx_tcga_vsd_scaled <- zScore(assay(pdx_tcga_vsd))

write.csv(pdx_tcga_vsd_scaled, "zScored_vst_pdx_")








