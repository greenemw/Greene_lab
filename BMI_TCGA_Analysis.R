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
colnames(tcga_bmi_cms)[c(1,4)] <- c("ID","cms")
tcga_bmi_cms$obese <- ifelse(tcga_bmi_cms$bmi >= 30, "Obese", "Not_Obese")

##Remove NA's in bmi data and add BMI category
dat <- tcga_bmi_cms[!is.na(dat$bmi),]
for (i in 1:length(dat$bmi)){
  if (dat$bmi[i] >=30){
    dat$bmiCat[i] <- "Obese"
  } else if(dat$bmi[i] < 30 & dat$bmi[i] >=25){
    dat$bmiCat[i] <- "Overweight"
  } else if(dat$bmi[i] < 25 & dat$bmi[i] >=18.5){
    dat$bmiCat[i] <- "Normal"
  } else (dat$bmiCat[i] <- NA)
}


write.csv(dat, "TCGA_BMI_with_Guinney_CMS_labels.csv")
##Explore BMI and CMS associations
attach(dat)
table(cms)/length(cms)
table(obese)/length(obese)
summary(bmi)

table(cms, bmiCat)

##Graphical Analysis
png("CMS and BMI Boxplots.png")
ggplot(data = dat, aes(cms,bmi)) + 
  geom_boxplot() + 
  labs(title = "BMI and CMS", x = "CMS", y = "BMI")
dev.off()

png("CMS and BMI Histograms.png")
ggplot(dat, aes(bmi)) + 
  geom_histogram(binwidth = 3,col = "black", fill = "white") +
  facet_wrap(~cms, ncol = 2) + 
  labs(title = "BMI and CMS Histograms", x = "BMI", y = "Frequency")
dev.off()

png("CMS and BMI eCDF.png")
ggplot(data = dat, aes(bmi, color = cms)) + 
  stat_ecdf(size = 1) + 
  labs(title = "Empirical Distribution of BMI", x = "BMI", y = "P(X)")
dev.off()

##Statistical Analysis
library(coin)

##Analyze with BMI as dichotomous
tblBmiCms <- table(obese, cms)
chisq.test(tblBmiCms)
independence_test(tblBmiCms, "quadratic")

##Analyze with BMI in 3 groups
tblBmiCatCms <- table(bmiCat,cms)
chisq.test(tblBmiCatCms)
independence_test(tblBmiCatCms, "quadratic")

##Analyze BMI as continuous
##Print summaries of each subtype
for (i in unique(cms)){
  print(paste("Summary for", i))
  print(summary(bmi[cms == i]))
}

independence_test(bmi~cms, alternative = "greater")
##Permutation test by hand
bmiFit <- lm(bmi ~ cms)
trueStat <- anova(bmiFit)[[4]][1]

R <- 5000
i <- 0
Fstat <- c()
while (i < R){
  permBmi <- base::sample(x = bmi,replace = F)
  permFit <- lm(permBmi ~ cms)
  Fstat[i] <- anova(permFit)[[4]][1]
  i <- i+1
}

Fn <- ecdf(Fstat)
1-Fn(trueStat)
png(filename = "ANOVA Permutation Test.png")
ggplot(data.frame(Fstat), aes(Fstat)) +
  geom_histogram(binwidth = 0.3, col="black", fill = "white") +
  geom_vline(xintercept = trueStat, linetype = "dashed", col = "black", size = 1) +
  labs(title = "Permutation ANOVA",
       subtitle = "Null Distribution of F",
        x = "Permuted F-statistic",
        y = "Count",
       caption = paste0("P-value for true F = ", round(1-Fn(trueStat),4)))
dev.off()

##Chi-sq permutation by hand
R <- 5000
i <- 0
Xstat <- c()
while(i < R) {
  permCms <- base::sample(cms, replace = F)
  permTbl <- table(obese,permCms)
  Xstat[i] <- chisq.test(permTbl)[[1]][[1]]
  i <- i+1
}

FnX <- ecdf(Xstat)
1-FnX(chisq$statistic)
hist(Xstat)


