##BMI Analysis of TCGA-COAD tissue samples
library(TCGAbiolinks)
library(SummarizedExperiment)
library(devtools)
library(tidyverse)

##Query the HT-Seq counts data from COAD project
query.RNAseq.COAD <- GDCquery(project = "TCGA-COAD", 
                              data.category = "Transcriptome Profiling",
                              data.type = "Gene Expression Quantification",
                              workflow.type = "HTSeq - Counts")
##Downloads the query above
GDCdownload(query.RNAseq.COAD)
##Puts all of the data together into a summarized experiment object. 
##Adds in the clinical data as well.
COADRNAseqData <- GDCprepare(query = query.RNAseq.COAD,
                             summarizedExperiment = TRUE)

##Put the TCGA gene counts data into a data frame for classification
TCGACOUNTS <- data.frame(assay(COADRNAseqData))
colnames(TCGACOUNTS) <- colnames(assay(COADRNAseqData))

##BMI data
tcgaBmi <- data.frame(sample = COADRNAseqData$patient, 
                      barcode = COADRNAseqData$barcode,
                      bmi = COADRNAseqData$bmi)


##Pulled CMS labels from Guinney et al.
cmsLabs <- read.table("cms_labels_public_all.txt",header = T)
tcgaCmsLabs <- subset(cmsLabs,
                      subset = dataset == "tcga", 
                      select = c(sample, CMS_final_network_plus_RFclassifier_in_nonconsensus_samples))

##Merge Guinney labels with BMI from TCGA
BMI_CMS <- merge(tcgaBmi,tcgaCmsLabs, by = "sample")
tcga_bmi_cms <- BMI_CMS[!duplicated(BMI_CMS$sample),]
colnames(tcga_bmi_cms)[4] <- "cms"
tcga_bmi_cms$obese <- ifelse(tcga_bmi_cms$bmi >= 30, "Obese", "Not Obese")

write.csv(tcga_bmi_cms, "TCGA_BMI_with_Guinney_CMS_labels.csv")

##Remove NA's in 
dat <- tcga_bmi_cms[!is.na(tcga_bmi_cms$bmi),-2]
##Explore BMI and CMS associations
attach(dat)
table(cms)/length(cms)
table(obese)/length(obese)
summary(bmi)

##Graphical Analysis
ggplot(data = dat, aes(cms,bmi)) + 
  geom_boxplot() + 
  labs(title = "BMI and CMS", x = "CMS", y = "BMI")
ggplot(dat, aes(bmi)) + 
  geom_histogram(binwidth = 3,col = "black", fill = "white") +
  facet_wrap(~cms, ncol = 2) + 
  labs(title = "BMI and CMS Histograms", x = "BMI", y = "Frequency")
ggplot(data = dat, aes(bmi, color = cms)) + 
  stat_ecdf(size = 1) + 
  labs(title = "Empirical Distribution of BMI", x = "BMI", y = "P(X)")

##Statistical Analysis
library(coin)
tblBmiCms <- table(obese, cms)
chisq.test(tblBmiCms)

independence_test(bmi~cms)




