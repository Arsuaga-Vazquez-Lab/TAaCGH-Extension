################################################################################
# This program downloads RNA-Seq, ABSOLUTE Copy number, Methylation 
# (Promoter region: TSS200, TSS1500, 5'UTR and 1stExon), Mutation,
# and miRNA data along with clinical information from TCGA-BRCA dataset.
# Gene and miRNA set can be modified by changing the list below. 
# Since methylation data is huge, it will require 16GB RAM. 
# Methylation section (boxed with double lines) can be skipped without any issue. 
# As we add more data, the sample size decreases due to limited data.
# TCGA_geneOfInterestTP holds gene expression (and clinical) data only.
# TCGA_geneOfInterestTP_CNV holds gene expression and CNV data.
# TCGA_geneOfInterestTP_CNV_miRNA holds gene expression, CNV, and miRNA data.
# After completing the program, TCGA_geneOfInterest holds all data. 
#
# This program should take about 10 min to fully load. Please run, one at a time,
# each section as indicated by a start comment "#== ..===" and ending comment
# "###...###". Program lines under "#=#...=#=" are optional. Recommended memory 
# RAM size of 16.0 Gb or more. Please email for any questions about this program.
# References: https://rdrr.io/bioc/TCGAbiolinks/
################################################################################
#======== Install libraries only once. Please restart RStudio after installing 
# the libraries ========
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("TCGAbiolinks", "SummarizedExperiment", "DESeq2", "sesame",
               "IlluminaHumanMethylation450kanno.ilmn12.hg19", "minfi", "multiMiR")
BiocManager::install(bioc_pkgs)

cran_pkgs <- c("MASS", "tidyverse", "survminer", "survival", "vioplot",
               "corrplot", "Hmisc", "ppcor", "glmnet", "DescTools",
               "randomForest", "vcd", "forcats", "shiny")
install.packages(cran_pkgs)
################################################################################

################################################################################
#======== Please load the libraries before every use ======== 
library(TCGAbiolinks)
library(tidyverse)
library(SummarizedExperiment)
library(survminer)
library(survival)
library(DESeq2)
library(vioplot)
library(corrplot)
library(sesame)
library(minfi)
library(MASS)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(Hmisc)
library(ppcor)
library(glmnet)
library(vcd)
library(forcats)
library(ggplot2)
library(dplyr)
################################################################################

################################################################################
#======== Set parameters to obtain LumA sample vector with desired genes========
subtype_data <- TCGAquery_subtype(tumor = "BRCA")

lumA <- list("PURPL", "TP53", "MYBBP1A", "BID", "CASP8", "MTOR", "ULK1",
             "PRKAA1", "RBM4", "ZBTB7A", "PIK3CA", "AKT1", "NFKB1", "PDK1",
             "RHEB", "CDC7", "NEK2", "RACGAP1") 
geneList <- c(lumA) #mainGene, p53RegGenes

miRNAList <- list("mir-137")
subtype <- c("LumA") #"LumA", "LumB", "Her2", "Basal", "Normal"
lumA_barcodes <- subtype_data$patient[subtype_data$BRCA_Subtype_PAM50 == subtype]
################################################################################

#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#======== Download TCGA-BRCA Methylation Data ======== 
#Please Skip this box (and next) if RAM is lower than 16GB
queryMethyl <- GDCquery(project = "TCGA-BRCA",
                        data.category = "DNA Methylation",
                        data.type = "Methylation Beta Value",
                        platform = "Illumina Human Methylation 450",
                        barcode = lumA_barcodes,
                        access = "open")
out_dir <- "~/Data/TCGA/GDC_BRCA_methyl"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = max(600, getOption("timeout")))  # give the network more time

GDCdownload(
  queryMethyl,
  method = "api",          # or "client" if you have the GDC client installed (see below)
  directory = out_dir,
  files.per.chunk = 10     # smaller tarballs -> fewer truncations
)

#GDCdownload(queryMethyl)
################################################################################

#Please skip this box (and prev.) if RAM is lower than 16GB
brcaMethyl <- GDCprepare(queryMethyl, summarizedExperiment = TRUE,directory = out_dir)
brcaMethylMatrix <- assay(brcaMethyl)
brcaMethylMatrix <- as.data.frame(brcaMethylMatrix)
write.csv(brcaMethylMatrix,"~/Data/TCGA/TCGA_methylation.csv")
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

################################################################################
#======== Download TCGA-BRCA RNA-Seq data ========
queryTP <- GDCquery(project = "TCGA-BRCA",
                    data.category = "Transcriptome Profiling",
                    experimental.strategy = "RNA-Seq",
                    workflow.type = "STAR - Counts",
                    barcode = lumA_barcodes,
                    access = "open")
GDCdownload(queryTP)
################################################################################

#======== Load the RNA-Seq data as summarized experiment and format data for 
# study ========
brcaTP <- GDCprepare(queryTP, summarizedExperiment = TRUE)
brcaTPMatrix <- as.data.frame(assay(brcaTP, "unstranded"))
brcaTPGeneMetadata <- as.data.frame(rowData(brcaTP))
clinicalData <- as.data.frame(colData(brcaTP))
################################################################################

#======== Normalize gene expression by apply the variance stabilizing 
# transformation to the raw RNA-Seq data ========
ddsTP <- DESeqDataSetFromMatrix(countData = brcaTPMatrix,
                                colData = clinicalData,
                                design = ~1)
vstTP <- vst(ddsTP, blind = FALSE)
brcaTPvst <- assay(vstTP)
brcaTPvst <- as.data.frame(brcaTPvst)
################################################################################

#======== Download ER, PR, and HER2 status ========
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Clinical",
                  data.type = "Clinical Supplement",
                  data.format = "BCR Biotab")
GDCdownload(query)
################################################################################

#======== Load the hormonal receptor data =========
clinical.all <- GDCprepare(query)
tcga_brca.clin <- as.data.frame(clinical.all$clinical_patient_brca)
tcga_brca.clin_original <- tcga_brca.clin
clinicalData_original <- clinicalData
tcga_brca.clin <- tcga_brca.clin[-c(1,2), ]# Remove the first two rows (metadata and column descriptions)
tcga_brca.clin <- tcga_brca.clin %>%
  dplyr::select(bcr_patient_barcode, er_status_by_ihc, pr_status_by_ihc, her2_status_by_ihc)
colnames(tcga_brca.clin) <- c("bcr_patient_barcode", "ER_Status", "PR_Status", "HER2_Status")
str(tcga_brca.clin)

# Load the clinical dataset from the summarized experiment
clinicalData <- as.data.frame(colData(brcaTP))
clinicalData$bcr_patient_barcode <- gsub("-[0-9][0-9][A-Z]$", "", clinicalData$bcr_patient_barcode)
clinicalData$bcr_patient_barcode <- trimws(toupper(clinicalData$bcr_patient_barcode))
tcga_brca.clin$bcr_patient_barcode <- trimws(toupper(tcga_brca.clin$bcr_patient_barcode))
clinicalData <- merge(clinicalData, tcga_brca.clin, by = "bcr_patient_barcode", all.x = TRUE)
################################################################################

#======== Download TCGA-BRCA Absolute Copy Number data ========
queryACN <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Copy Number Variation",
                     data.type = "Gene Level Copy Number",
                     workflow.type = "ABSOLUTE LiftOver",
                     #barcode = lumA_barcodes,
                     access = "open")
GDCdownload(queryACN)
################################################################################

#======== Load the Copy Number data as summarized experiment and format data for 
# study ========
brcaACN <- GDCprepare(queryACN, summarizedExperiment = TRUE)
brcaACNMatrix <- assay(brcaACN)
brcaACNMatrix <- as.data.frame(brcaACNMatrix)
brcaACNMetadata <- as.data.frame(rowData(brcaACN))
################################################################################

#======== Download TCGA-BRCA Mutation data ========
queryMut <- GDCquery(project = "TCGA-BRCA",
                     data.category = "Simple Nucleotide Variation",
                     data.type = "Masked Somatic Mutation",
                     workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
                     barcode = lumA_barcodes,
                     access = "open")
GDCdownload(queryMut)

#Load the Mutation data as summarized experiment and format data for study
brcaMut <- GDCprepare(queryMut, summarizedExperiment = TRUE)
################################################################################
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#======== Download TCGA-BRCA miRNA data ========
querymiRNA <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "miRNA Expression Quantification",
                       workflow.type = "BCGSC miRNA Profiling",
                       barcode = lumA_barcodes,
                       sample.type = c("Primary Tumor", "Solid Tissue Normal"))
GDCdownload(querymiRNA)
################################################################################

#======== Load the miRNA data as summarized experiment and format data for 
# study ========
brcamiRNA<- GDCprepare(querymiRNA)
miRNA_IDs <- brcamiRNA[, "miRNA_ID"]
brcamiRNA_RPM <- brcamiRNA[, grepl("reads_per_million", colnames(brcamiRNA))]
brcamiRNA_RPM<- cbind(miRNA_IDs, brcamiRNA_RPM)
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
################################################################################

#Extracted data from TCGA above ^

#Now use extracted data to create a summarized table below v

#======== Create overall survival information using vital status, days to last 
# followup, and days to death ========
clinicalData$overallSurvival <- ifelse(clinicalData$paper_vital_status == "Alive",
                                       clinicalData$paper_days_to_last_followup,
                                       clinicalData$paper_days_to_death) #Alive is censored case
clinicalData<- clinicalData[!is.na(clinicalData$overallSurvival),] #Remove missing time
clinicalData$deceased <- ifelse(clinicalData$paper_vital_status== "Alive", 0, 1) #0=censored, 1=deaths
participantCode <- sub(".*-(....)$", "\\1", clinicalData$patient)
duplicates <- participantCode[duplicated(participantCode)]
#print(duplicates)
uniqueParticipantCode <- make.unique(participantCode)
clinicalData$participantCode <- uniqueParticipantCode
TCGA_geneOfInterest <- clinicalData[!grepl("\\.[1-9]$", clinicalData$participantCode), ]#debug remove duplicates
################################################################################

#======== Function Indexed Reference Lists for Genes and miRNAs ========
assign_genes <- function(...) {
  genes <- c(...)
  gene_list <- list()
  for (i in seq_along(genes)) {
    gene_list[[paste0("geneName", i)]] <- genes[[i]]
    gene_list[[paste0("geneNameExp", i)]] <- paste(genes[[i]], "Exp", sep = "_")
  }
  return(gene_list)
}
assign_miRNA <- function(...){
  miRNAs <- c(...)
  miRNA_list <- list()
  for (i in seq_along(miRNAs)){
    miRNA_list[[paste0("miRNA",i)]] <- miRNAs[[i]]
  }
  return(miRNA_list)
}
geneList <- assign_genes(geneList)
miRNAList <- assign_miRNA(miRNAList)
################################################################################

#======== Function to add gene expression to main dataset =========
column_names <- colnames(brcaTPvst)
participantCodesTP <- substr(column_names, 9, 12)
colnames(brcaTPvst) <- make.unique(participantCodesTP)
cols_to_remove <- grep("\\.[0-9]$", colnames(brcaTPvst))
cols_with_digit_suffix <- grep("\\.[0-9]$", colnames(brcaTPvst))
num_cols_with_digit_suffix <- length(cols_with_digit_suffix)
print(num_cols_with_digit_suffix)
brcaTPvst <- brcaTPvst[, -cols_to_remove]

addGeneExp <- function(geneName) {
  TP_geneOfInterest <- brcaTPvst %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    gather(key = "case_id", value = "counts", -gene_id) %>%
    left_join(., brcaTPGeneMetadata, by = "gene_id") %>%
    filter(gene_name == geneName)
  tempVars <- c("case_id", "counts")
  TP_geneOfInterest <- TP_geneOfInterest[tempVars]
  colnames(TP_geneOfInterest)[colnames(TP_geneOfInterest) == "case_id"] <- "participantCode"
  geneNameExp <- paste(geneName, "Exp", sep = "_")
  colnames(TP_geneOfInterest)[2] <- geneNameExp
  if (nrow(TP_geneOfInterest) > 0) {
    # Perform the merge only if TP_geneOfInterest is not empty
    TCGA_geneOfInterest <- merge(TP_geneOfInterest, TCGA_geneOfInterest, by = "participantCode")
  }
  return(TCGA_geneOfInterest)
}
################################################################################

#======== Apply the function to get gene expression onto main dataset ======== 
for (i in seq_along(geneList)) {
  geneKey <- paste0("geneName", i)
  if (geneKey %in% names(geneList)) {
    geneName <- geneList[[geneKey]]
    tryCatch({
      TCGA_geneOfInterest <- addGeneExp(geneName)
    }, error = function(e) {
      message("Skipping ", geneName, " due to error: ", e$message)
    })
  }
}

#TCGA_geneOfInterestTP <- TCGA_geneOfInterest #Creates a copy of table with GE for debugging

################################################################################

#======== Function to add Absolute Copy Number to main dataset =========
column_names <- colnames(brcaACNMatrix)
participantCodesACN <- substr(column_names, 9, 12)
colnames(brcaACNMatrix) <- make.unique(participantCodesACN)
cols_to_remove <- grep("\\.[0-9]$", colnames(brcaACNMatrix))
num_cols_to_remove <- length(cols_to_remove)
print(num_cols_to_remove)
brcaACNMatrix <- brcaACNMatrix[, -cols_to_remove]

getACN <- function(geneName) {
  ACN_geneOfInterest <- brcaACNMatrix %>%
    as.data.frame() %>%
    rownames_to_column(var = "gene_id") %>%
    gather(key = "case_id", value = "counts", -gene_id) %>%
    left_join(., brcaACNMetadata, by = "gene_id") %>%
    filter(gene_name == geneName)
  geneNameACN <- paste(geneName, "ACN", sep = "_")
  colnames(ACN_geneOfInterest) <- c("gene_id", "case_id", geneNameACN, "gene_name")
  myvars <- c("case_id", geneNameACN)
  ACN_geneOfInterest <- ACN_geneOfInterest[myvars]
  colnames(ACN_geneOfInterest)[colnames(ACN_geneOfInterest) == "case_id"] <- "participantCode"
  if (nrow(ACN_geneOfInterest) > 0) {
    TCGA_geneOfInterest <- merge(ACN_geneOfInterest, TCGA_geneOfInterest, by = "participantCode")
  }
  return(TCGA_geneOfInterest)
}
################################################################################

#======== Apply the function to get copy number onto main dataset ======== 
for (i in seq_along(geneList)) {
  geneKey <- paste0("geneName", i)
  if (geneKey %in% names(geneList)) {
    geneName <- geneList[[geneKey]]
    TCGA_geneOfInterest <-getACN(geneName)
  }
}

#TCGA_geneOfInterestTP_CNV <- TCGA_geneOfInterest #Creates a copy of table with GE+CN for debugging
#TCGA_geneOfInterest <- TCGA_geneOfInterestTP_CNV 

################################################################################
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=

#======== Function to add Micro RNA to main dataset =========
add_miRNA_exp <- function(miRNA_namess) {
  # Check if the miRNA is in the dataset
  miRNA_row <- brcamiRNA_RPM[brcamiRNA_RPM$miRNA_ID == miRNA_namess, ]
  
  if (nrow(miRNA_row) == 0) {
    message("miRNA ", miRNA_names, " not found in dataset.")
    return(TCGA_geneOfInterest)
  }
  miRNA_exp <- t(miRNA_row[, -1])  # drop the first column (miRNA_ID)
  miRNA_exp <- as.data.frame(miRNA_exp)
  sample_names <- colnames(brcamiRNA_RPM)[-1]  # exclude 'miRNA_ID'
  rownames(miRNA_exp) <- sample_names
  miRNA_exp$participantCode <- sapply(sample_names, function(x) strsplit(x, "-")[[1]][3])
  colnames(miRNA_exp)[1] <- paste0(miRNA_names, "_miRNAExp")
  
  # Merge into TCGA_geneOfInterest
  merged_data <- merge(TCGA_geneOfInterest, miRNA_exp[, c("participantCode", paste0(miRNA_names, "_miRNAExp"))],
                       by = "participantCode", all.x = TRUE)
}

################################################################################

#======== Apply the function to get miRNA onto main dataset ======== 
for (i in seq_along(miRNAList)) {
  miRNAKey <- paste0("miRNA", i)
  
  if (miRNAKey %in% names(miRNAList)) {
    miRNA_names <- miRNAList[[miRNAKey]]
    
    tryCatch({
      TCGA_geneOfInterest <- add_miRNA_exp(miRNA_names)
    }, error = function(e) {
      message("Skipping ", miRNA_names, " due to error: ", e$message)
    })
  }
}
#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=#=
################################################################################

#======== Function to add Mutation to main dataset =========
brcaMut$participantCode <- substr(brcaMut$Tumor_Sample_Barcode, 9, 12)
getMutVector <- function(geneName) {
  Mut_geneOfInterest <- brcaMut %>% 
    filter(Hugo_Symbol == geneName)
  MutVector <- Mut_geneOfInterest$participantCode
  return(MutVector)
}

#======== Apply the function to get Mutation onto main dataset ======== 
for (i in seq_along(geneList)) {
  geneKey <- paste0("geneName", i)
  if (geneKey %in% names(geneList)) {
    geneName <- geneList[[geneKey]]
    mutVector <- getMutVector(geneName)
    if (length(mutVector) > 0) {
      geneNameMutStatus <- paste0(geneName, "_Mut_Status")
      geneNameMut <- paste0(geneName, "_Mut")
      geneNameWT <- paste0(geneName, "_WT")
      TCGA_geneOfInterest[[geneNameMutStatus]] <- ifelse(
        TCGA_geneOfInterest$participantCode %in% mutVector, geneNameMut, geneNameWT)
    }
  }
}
TCGA_geneOfInterest <- TCGA_geneOfInterest[!duplicated(TCGA_geneOfInterest$participantCode), ]
#Ensure each patient has Mut status attached
TCGA_geneOfInterest <- TCGA_geneOfInterest[, apply(TCGA_geneOfInterest, 2, 
                                                   function(x) length(unique(x)) > 1)] 
################################################################################
#This completes the program. The dataset is in "TCGA_geneOfInterest"
################################################################################
TCGA_flat <- tibble::as_tibble(lapply(TCGA_geneOfInterest, function(x) if (is.list(x)) sapply(x, toString) else x))
#write.csv(TCGA_flat, "~/Data/TCGA/TCGA_LumA.csv", row.names = FALSE)
################################################################################

