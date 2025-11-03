################################################################################
df <- read.csv("~/Data/SCAN-B/SCAN-B_GSE96058.csv", stringsAsFactors = FALSE)
TCGA_geneOfInterest <- df
################################################################################

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

################################################################################
########### Survival p values and adjusted p values (BH) using Exp. ############

# === Helper: per-group MaxStat + log-rank p (LumA only, no truncation) ===
process_expression_group <- function(df, group_label, group_type) {
  df <- df %>%
    dplyr::filter(
      pam50_subtype == "LumA",
      !is.na(PURPL_Exp),
      !is.na(overallSurvival), !is.na(deceased)
    ) %>%
    dplyr::mutate(
      overallSurvival = suppressWarnings(as.numeric(as.character(overallSurvival))),
      deceased        = suppressWarnings(as.numeric(as.character(deceased)))
    ) %>%
    dplyr::filter(!is.na(overallSurvival), !is.na(deceased))
  
  if (nrow(df) < 2 || length(unique(df$deceased)) < 2) {
    return(data.frame(
      Group_Type        = group_type,
      Group             = group_label,
      n                 = nrow(df),
      PURPL_Cutoff      = NA_real_,
      MaxStat_statistic = NA_real_,
      PURPL_p_value     = NA_real_
    ))
  }
  
  cutpoint <- tryCatch(
    surv_cutpoint(df, time = "overallSurvival", event = "deceased", variables = "PURPL_Exp"),
    error = function(e) NULL
  )
  if (is.null(cutpoint) || is.infinite(cutpoint$cutpoint[[1]]) || is.na(cutpoint$cutpoint[[1]])) {
    return(data.frame(
      Group_Type        = group_type,
      Group             = group_label,
      n                 = nrow(df),
      PURPL_Cutoff      = NA_real_,
      MaxStat_statistic = NA_real_,
      PURPL_p_value     = NA_real_
    ))
  }
  
  cutoff_val   <- cutpoint$cutpoint[[1]]
  maxstat_stat <- summary(cutpoint)$statistic[[1]]
  
  df <- df %>% dplyr::mutate(PURPL_HighLow = ifelse(PURPL_Exp <= cutoff_val, "Low", "High"))
  
  pval <- tryCatch({
    surv_obj <- Surv(df$overallSurvival, df$deceased)
    fit      <- survdiff(surv_obj ~ PURPL_HighLow, data = df)
    pchisq(fit$chisq, df = 1, lower.tail = FALSE)
  }, error = function(e) NA_real_)
  
  data.frame(
    Group_Type        = group_type,
    Group             = group_label,
    n                 = nrow(df),
    PURPL_Cutoff      = round(cutoff_val, 3),
    MaxStat_statistic = round(maxstat_stat, 4),
    PURPL_p_value     = pval
  )
}

# === Build combos (PIK3CA+TP53, PIK3CA+AKT1) + singles (PIK3CA, TP53, AKT1) ===
combo_gene_pairs <- list(c("PIK3CA","TP53"), c("PIK3CA","AKT1"))
combo_results <- list()
for (pair in combo_gene_pairs) {
  gene1 <- pair[1]; gene2 <- pair[2]
  col1 <- paste0(gene1, "_Mut_Status")
  col2 <- paste0(gene2, "_Mut_Status")
  if (!all(c(col1, col2) %in% colnames(df))) next
  for (s1 in c(paste0(gene1,"_WT"), paste0(gene1,"_Mut"))) {
    for (s2 in c(paste0(gene2,"_WT"), paste0(gene2,"_Mut"))) {
      df_sub <- df %>% dplyr::filter(!!sym(col1) == s1, !!sym(col2) == s2)
      label  <- paste(s1, "+", s2)
      res    <- process_expression_group(df_sub, label, paste0(gene1, " + ", gene2, " Combo"))
      if (!is.null(res)) combo_results[[length(combo_results)+1]] <- res
    }
  }
}

single_results <- list()
for (g in c("PIK3CA","TP53","AKT1")) {
  col <- paste0(g, "_Mut_Status"); if (!col %in% colnames(df)) next
  for (status in c(paste0(g,"_WT"), paste0(g,"_Mut"))) {
    df_sub <- df %>% dplyr::filter(!!sym(col) == status)
    res    <- process_expression_group(df_sub, status, "Single Mutation Status")
    if (!is.null(res)) single_results[[length(single_results)+1]] <- res
  }
}

global_result <- process_expression_group(df, "All Samples", "Overall Expression Effect")

# === Keep only these groups for BH ===
keep_groups <- c(
  "PIK3CA_WT + TP53_WT",
  "PIK3CA_WT + TP53_Mut",
  "PIK3CA_Mut + TP53_WT",
  "PIK3CA_Mut + TP53_Mut",
  "PIK3CA_WT + AKT1_WT",
  "PIK3CA_Mut + AKT1_WT",
  "PIK3CA_WT",
  "PIK3CA_Mut",
  "TP53_WT",
  "TP53_Mut",
  "AKT1_WT"
)

all_strata <- dplyr::bind_rows(combo_results, single_results)

expr_tbl_core <- all_strata %>%
  dplyr::filter(Group %in% keep_groups) %>%
  dplyr::mutate(.ord = match(Group, keep_groups)) %>%
  dplyr::arrange(.ord) %>%
  dplyr::select(-.ord)

expr_tbl_core$Adjusted_p <- NA_real_
mask <- !is.na(expr_tbl_core$PURPL_p_value)
if (any(mask)) {
  expr_tbl_core$Adjusted_p[mask] <-
    p.adjust(as.numeric(expr_tbl_core$PURPL_p_value[mask]), method = "BH")
}

expr_tbl <- expr_tbl_core
if (is.data.frame(global_result) && nrow(global_result) == 1) {
  global_row <- global_result
  global_row$Adjusted_p <- NA_real_  # not part of BH family
  expr_tbl <- dplyr::bind_rows(expr_tbl, global_row)
}
print(expr_tbl)
# write.csv(expr_tbl, "~/Data/SCAN-B/SCAN-B_pval_exp.csv", row.names = FALSE)
################################################################################





################################################################################
########### Expression Survival using MaxStat adjusted p-value (BH) ############

# ======= Select wanted group by uncommenting =======
df_km <- df %>%
  dplyr::filter(
    pam50_subtype == "LumA",
    #TP53_Mut_Status   == "TP53_WT",
    #PIK3CA_Mut_Status == "PIK3CA_WT",
    #AKT1_Mut_Status   == "AKT1_WT",
    !is.na(PURPL_Exp),
    !is.na(overallSurvival), !is.na(deceased)
  ) %>%
  dplyr::mutate(
    overallSurvival = as.numeric(overallSurvival),
    deceased        = as.numeric(deceased)
  ) %>%
  dplyr::filter(!is.na(overallSurvival), !is.na(deceased))

cutpoint <- tryCatch(
  surv_cutpoint(df_km, time = "overallSurvival", event = "deceased", variables = "PURPL_Exp"),
  error = function(e) NULL
)

if (is.null(cutpoint)) {
  plot.new(); title("MaxStat could not determine a cutoff for PURPL_Exp")
} else {
  df_km <- df_km %>%
    dplyr::mutate(PURPL_HighLow = ifelse(PURPL_Exp <= cutpoint$cutpoint[[1]], "Low", "High"))
  
  surv_obj <- Surv(df_km$overallSurvival, df_km$deceased)
  fit      <- survfit(surv_obj ~ PURPL_HighLow, data = df_km)
  raw_p    <- as.numeric(surv_pvalue(fit, data = df_km)$pval)
  
  mut_cols <- c("PIK3CA_Mut_Status","TP53_Mut_Status","AKT1_Mut_Status")
  mut_cols <- mut_cols[mut_cols %in% colnames(df_km)]
  
  fixed <- lapply(mut_cols, function(cl) {
    u <- unique(df_km[[cl]]); u <- u[!is.na(u)]
    if (length(u) == 1) as.character(u) else NA_character_
  })
  names(fixed) <- mut_cols
  fixed <- Filter(Negate(is.na), fixed)
  
  lookup_rows <- NULL
  if (length(fixed) >= 2) {
    genes <- sub("_Mut_Status$","", names(fixed))
    g1 <- genes[1]; g2 <- genes[2]
    s1 <- fixed[[1]]; s2 <- fixed[[2]]
    group_candidates <- c(paste0(s1, " + ", s2), paste0(s2, " + ", s1))
    type_candidates  <- c(paste0(g1, " + ", g2, " Combo"), paste0(g2, " + ", g1, " Combo"))
    lookup_rows <- subset(expr_tbl, Group %in% group_candidates & Group_Type %in% type_candidates)
  } else if (length(fixed) == 1) {
    s1 <- unname(fixed[[1]])
    lookup_rows <- subset(expr_tbl, Group_Type == "Single Mutation Status" & Group == s1)
  } else {
    lookup_rows <- subset(expr_tbl, Group_Type == "Overall Expression Effect" & Group == "All Samples")
  }
  
  disp_p <- raw_p
  if (is.data.frame(lookup_rows) && nrow(lookup_rows) >= 1) {
    p_bh  <- suppressWarnings(as.numeric(lookup_rows$Adjusted_p[1]))
    p_raw <- suppressWarnings(as.numeric(lookup_rows$PURPL_p_value[1]))
    if (!is.na(p_bh))      disp_p <- p_bh
    else if (!is.na(p_raw)) disp_p <- p_raw
  }
  pval_text <- paste0("p = ", signif(disp_p, 4))
  
  g <- ggsurvplot(
    fit,
    data = df_km,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.title = NULL,
    conf.int = FALSE,
    legend = c(0.85, 0.45),
    legend.title = "PURPL gene exp.",
    legend.labs = c("High", "Low"),
    legend.title.font = 2,
    legend.text.font = 2,
    xlab = "Time (days)",
    ylab = "Survival probability",
    font.x = c(30, "bold"),
    font.y = c(30, "bold"),
    font.tickslab = 30,
    font.legend = c(20,"bold"),
    tables.theme = theme_classic(base_size = 25),
    ggtheme = theme_classic() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color = "black", fill = "transparent", size = 0.5)
      ),
    pval = pval_text,
    pval.size = 8
  )
  
  # keep colored tags; remove redundant titles on risk table
  g$table <- g$table + labs(x = NULL, y = NULL, title = NULL) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  print(g)
}
################################################################################


