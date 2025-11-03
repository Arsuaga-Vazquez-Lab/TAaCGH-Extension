################################################################################
#======== Please load the needed libraries before every use ======== 
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

#Choose what dataset to compute.
TCGA_geneOfInterest<- read.csv("~/Data/TCGA/TCGA_LumA.csv") #TCGA
#TCGA_geneOfInterest<- read.csv("~/Data/SCAN-B/SCAN-B_GSE96058.csv") #SCANB

df<-TCGA_geneOfInterest

# == HELPER FUNCTION (carefully select filter for df based on dataset chosen) ==
process_expression_group <- function(df, group_label, group_type,
                                     pmethod = c("condMC")) { #Pick method "HL","Lau94","Lau92","condMC","exactGauss"
  pmethod <- match.arg(pmethod)
  
  df <- df %>%
    filter(
      !is.na(PURPL_Exp),
      !is.na(PURPL_ACN), !PURPL_ACN %in% c(0, 1, 7), ###### Uncomment if TCGA dataset, comment if SCAN-B
      !is.na(overallSurvival), !is.na(deceased),
      #pam50_subtype == "LumA" ##### Uncomment if SCAN-B dataset, comment if TCGA
    ) %>%
    mutate(
      overallSurvival = suppressWarnings(as.numeric(as.character(overallSurvival))),
      deceased        = suppressWarnings(as.numeric(as.character(deceased)))
    ) %>%
    filter(!is.na(overallSurvival), !is.na(deceased))
  
  # Need at least 2 events and 2 groups possible
  if (nrow(df) < 10 || length(unique(df$deceased)) < 2) {
    return(data.frame(
      Group_Type = group_type, Group = group_label, n = nrow(df),
      PURPL_Cutoff = NA, MaxStat_statistic = NA, PURPL_p_value = NA,
      Pmethod = pmethod, stringsAsFactors = FALSE
    ))
  }
  
  # Maxstat test: finds optimal cutpoint AND computes p-value corrected for the search
  mt <- tryCatch(
    maxstat::maxstat.test(
      Surv(overallSurvival, deceased) ~ PURPL_Exp,
      data    = df,
      smethod = "LogRank",
      pmethod = pmethod
    ),
    error = function(e) NULL
  )
  
  if (is.null(mt)) {
    return(data.frame(
      Group_Type = group_type, Group = group_label, n = nrow(df),
      PURPL_Cutoff = NA, MaxStat_statistic = NA, PURPL_p_value = NA,
      Pmethod = pmethod, stringsAsFactors = FALSE
    ))
  }
  
  cutoff_val   <- unname(mt$estimate)      # estimated cutpoint
  M_val        <- unname(mt$statistic)     # maximally selected (standardized) statistic
  pval_adj_cut <- unname(mt$p.value)       # p-value adjusted for scanning cutpoints
  
  data.frame(
    Group_Type = group_type,
    Group = group_label,
    n = nrow(df),
    PURPL_Cutoff = round(cutoff_val, 3),
    MaxStat_statistic = round(M_val, 4),
    PURPL_p_value = pval_adj_cut,   # already corrected for cutpoint search
    Pmethod = pmethod,
    stringsAsFactors = FALSE
  )
}

# ==== PARTS 1–4 (unchanged structure), now using the revised helper ====
combo_gene_pairs <- list(c("PIK3CA", "TP53"), c("PIK3CA"))
combo_results <- list()

for (pair in combo_gene_pairs) {
  gene1 <- pair[1]; gene2 <- pair[2]
  col1 <- paste0(gene1, "_Mut_Status"); col2 <- paste0(gene2, "_Mut_Status")
  if (!all(c(col1, col2) %in% colnames(TCGA_geneOfInterest))) next
  
  for (status1 in c(paste0(gene1, "_WT"), paste0(gene1, "_Mut"))) {
    for (status2 in c(paste0(gene2, "_WT"), paste0(gene2, "_Mut"))) {
      df_sub <- TCGA_geneOfInterest %>% filter(!!rlang::sym(col1) == status1,
                                               !!rlang::sym(col2) == status2)
      label <- paste(status1, "+", status2)
      res <- process_expression_group(df_sub, label, paste0(gene1, " + ", gene2, " Combo"))
      if (!is.null(res)) combo_results[[length(combo_results) + 1]] <- res
    }
  }
}

single_results <- list()
mutation_genes1 <- unique(c(mutation_genes, sub("_.*", "", driver_statuses)))
statuses_to_test <- unlist(lapply(mutation_genes1, function(g) c(paste0(g, "_WT"), paste0(g, "_Mut"))))

for (status in statuses_to_test) {
  gene <- sub("_.*", "", status)
  col <- paste0(gene, "_Mut_Status")
  if (!col %in% colnames(TCGA_geneOfInterest)) next
  df_sub <- TCGA_geneOfInterest %>% filter(!!rlang::sym(col) == status)
  res <- process_expression_group(df_sub, status, "Single Mutation Status")
  if (!is.null(res)) single_results[[length(single_results) + 1]] <- res
}

global_result <- process_expression_group(
  TCGA_geneOfInterest, group_label = "All Samples", group_type = "Overall Expression Effect"
)

# ==== COMBINE ====
expr_final <- dplyr::bind_rows(combo_results, single_results)

# Append global result (PURPL background)
if (!is.null(global_result)) expr_final <- dplyr::bind_rows(expr_final, global_result)

print(expr_final)

#write.csv(expr_final,"~/Data/TCGA/TCGA_LumA_maxstat_pvalues.csv") #TCGA
#write.csv(expr_final,"~/Data/SCAN-B/TCGA_LumA_maxstat_pvalues.csv") #SCAN-B




