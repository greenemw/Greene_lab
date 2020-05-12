# RNASeq_Analysis
Analysis pipeline for RNA seq data from the Greene Lab

The pipeline currently will run three different analyses based in three different packages: DESeq2, CMScaller, and TCGAbiolinks

The pipeline will take in a RNA seq count file and and a phenotype table for the experiment and run DESeq2. 
It will output a PCA plot of the experimental covariates and contrasts specified by the user.

The CMScaller portion of the pipeline will output a table of assigned subtypes and a gene set analysis.

The TCGAbiolinks portion will run a DESeq2 analysis against the normal tissue from the TCGA-COAD project.
Similar DESeq contrasts can be performed with this portion of the pipeline.
