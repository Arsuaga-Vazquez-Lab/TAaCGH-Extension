################################################################################
df <- read.csv("~/Data/TCGA/TCGA_LumA.csv",stringsAsFactors = FALSE)
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
# Survival p values and adjusted p values (BH) using Exp. Truncation 2500 days #

#Truncating at 2500 for KM survival TCGA expression:
mutation_genes <- c("AKT1")
driver_statuses <- c("PIK3CA_WT", "PIK3CA_Mut", "TP53_WT", "TP53_Mut")
combine_statuses <- FALSE  # set to TRUE if needed
TAU <- 2500  # days

# ---- Helper to cap follow-up at TAU ----
cap_followup <- function(df, time_col = "overallSurvival", event_col = "deceased", tau = TAU) {
  time  <- suppressWarnings(as.numeric(as.character(df[[time_col]])))
  event <- suppressWarnings(as.numeric(as.character(df[[event_col]])))
  ok <- !is.na(time) & !is.na(event)
  time_cap  <- time
  event_cap <- event
  time_cap[ok]  <- pmin(time[ok], tau)
  # anyone beyond tau becomes censored at tau
  event_cap[ok & time[ok] > tau] <- 0
  df[[time_col]]  <- time_cap
  df[[event_col]] <- event_cap
  df
}

# ==== HELPER FUNCTION (with truncation at TAU) ====
process_expression_group <- function(df, group_label, group_type, tau = TAU) {
  df <- df %>%
    dplyr::filter(
      !is.na(PURPL_Exp),
      !is.na(PURPL_ACN), !PURPL_ACN %in% c(0, 1, 7),
      !is.na(overallSurvival), !is.na(deceased)
    ) %>%
    dplyr::mutate(
      overallSurvival = suppressWarnings(as.numeric(as.character(overallSurvival))),
      deceased        = suppressWarnings(as.numeric(as.character(deceased)))
    ) %>%
    dplyr::filter(!is.na(overallSurvival), !is.na(deceased)) %>%
    cap_followup(time_col = "overallSurvival", event_col = "deceased", tau = tau)
  
  # after truncation, ensure there are events on or before tau
  if (nrow(df) < 2 || length(unique(df$deceased)) < 2 || sum(df$deceased == 1, na.rm = TRUE) == 0) {
    return(data.frame(
      Group_Type        = group_type,
      Group             = group_label,
      n                 = nrow(df),
      PURPL_Cutoff      = NA_real_,
      MaxStat_statistic = NA_real_,
      PURPL_p_value     = NA_real_
    ))
  }
  
  # cutpoint uses truncated time/event
  cutpoint <- tryCatch({
    surv_cutpoint(df, time = "overallSurvival", event = "deceased", variables = "PURPL_Exp")
  }, error = function(e) NULL)
  
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
    surv_obj <- survival::Surv(df$overallSurvival, df$deceased)
    fit <- survival::survdiff(surv_obj ~ PURPL_HighLow, data = df)
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


# ==== PART 1: DRIVER + MUTATION COMBINATIONS ====
combo_gene_pairs <- list(
  c("PIK3CA", "TP53"),
  c("PIK3CA", "AKT1")
)

combo_results <- list()

for (pair in combo_gene_pairs) {
  gene1 <- pair[1]
  gene2 <- pair[2]
  
  col1 <- paste0(gene1, "_Mut_Status")
  col2 <- paste0(gene2, "_Mut_Status")
  
  if (!all(c(col1, col2) %in% colnames(TCGA_geneOfInterest))) next
  
  for (status1 in c(paste0(gene1, "_WT"), paste0(gene1, "_Mut"))) {
    for (status2 in c(paste0(gene2, "_WT"), paste0(gene2, "_Mut"))) {
      
      df_sub <- TCGA_geneOfInterest %>%
        filter(!!sym(col1) == status1, !!sym(col2) == status2)
      
      label <- paste(status1, "+", status2)
      
      result <- process_expression_group(df_sub, label, paste0(gene1, " + ", gene2, " Combo"))
      if (!is.null(result)) combo_results[[length(combo_results) + 1]] <- result
    }
  }
}


# ==== PART 2: SINGLE MUTATION STATUS ====
single_results <- list()
mutation_genes1 <- unique(c(mutation_genes, sub("_.*", "", driver_statuses)))
statuses_to_test <- unlist(lapply(mutation_genes1, function(g) c(paste0(g, "_WT"), paste0(g, "_Mut"))))

for (status in statuses_to_test) {
  gene <- sub("_.*", "", status)
  col <- paste0(gene, "_Mut_Status")
  if (!col %in% colnames(TCGA_geneOfInterest)) next
  
  df_sub <- TCGA_geneOfInterest %>% filter(!!sym(col) == status)
  result <- process_expression_group(df_sub, status, "Single Mutation Status")
  if (!is.null(result)) single_results[[length(single_results) + 1]] <- result
}

# ==== PART 3: GLOBAL EXPRESSION SURVIVAL ====
global_result <- process_expression_group(
  TCGA_geneOfInterest,
  group_label = "All Samples",
  group_type = "Overall Expression Effect"
)

# ==== PART 4: CUSTOM TRIPLE-MUTATION GROUPS ====
triple_combo_results <- list()

# Define the three relevant mutation status columns
status_cols <- c("PIK3CA_Mut_Status", "AKT1_Mut_Status", "MTOR_Mut_Status")

# Check if all required columns exist in data
if (all(status_cols %in% colnames(TCGA_geneOfInterest))) {
  # Define your two desired triple groups
  triple_groups <- list(
    list(
      label = "PIK3CA_WT + AKT1_WT + MTOR_WT",
      filters = c("PIK3CA_WT", "AKT1_WT", "MTOR_WT")
    ),
    list(
      label = "PIK3CA_WT + AKT1_WT + MTOR_Mut",
      filters = c("PIK3CA_WT", "AKT1_WT", "MTOR_Mut")
    )
  )
  
  for (group_def in triple_groups) {
    filter_vals <- group_def$filters
    group_label <- group_def$label
    
    df_sub <- TCGA_geneOfInterest %>%
      filter(
        PIK3CA_Mut_Status == filter_vals[1],
        AKT1_Mut_Status == filter_vals[2],
        MTOR_Mut_Status == filter_vals[3]
      )
    
    result <- process_expression_group(df_sub, group_label, "PIK3CA + AKT1 + MTOR Combo")
    if (!is.null(result)) triple_combo_results[[length(triple_combo_results) + 1]] <- result
  }
}

# ==== COMBINE AND ADJUST P-VALUES ====
expr_final <- bind_rows(combo_results, single_results, triple_combo_results)

# Adjust across all groups
expr_final <- expr_final %>%
  mutate(Adjusted_p = p.adjust(PURPL_p_value, method = "BH"))


# Append global result (what PURPL expression is)
if (!is.null(global_result)) {
  expr_final <- bind_rows(expr_final, global_result)
}

# ==== OUTPUT ====
print(expr_final)
# ==== KEEP ONLY TABLE ROWS & DO BH JUST ON THEM ====

# The exact groups we want in the final table (order preserved below)
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
  "AKT1_WT",
  "All Samples"
)

# Keep only those rows (and keep the table order)
expr_tbl <- expr_final %>%
  dplyr::filter(Group %in% keep_groups) %>%
  dplyr::mutate(.order = match(Group, keep_groups)) %>%
  dplyr::arrange(.order) %>%
  dplyr::select(-.order)

# Compute BH only across rows that should be adjusted
# (exclude "All Samples" and NA p-values)
idx_adjust <- with(expr_tbl, Group != "All Samples" & !is.na(PURPL_p_value))

expr_tbl$Adjusted_p <- NA_real_
if (any(idx_adjust)) {
  expr_tbl$Adjusted_p[idx_adjust] <-
    p.adjust(as.numeric(expr_tbl$PURPL_p_value[idx_adjust]), method = "BH")
}

#set p-values to 4 digits
expr_tbl <- expr_tbl %>%
  dplyr::mutate(
    PURPL_Cutoff = round(PURPL_Cutoff, 3),
    MaxStat_statistic = round(MaxStat_statistic, 4),
    PURPL_p_value = as.numeric(PURPL_p_value),   # ensure numeric
    Adjusted_p = ifelse(is.na(Adjusted_p), NA, Adjusted_p)
  )

# ==== OUTPUT ====
print(expr_tbl)
#write.csv(expr_tbl,"~/Data/TCGA/TCGA_pval_exp.csv")
################################################################################





################################################################################
########### Survival p values and adjusted p values (BH) using ACN #############

# ==== PARAMETERS ====
mutation_genes <- c("AKT1")
driver_statuses <- c("PIK3CA_WT", "PIK3CA_Mut", "TP53_WT", "TP53_Mut")
combine_statuses <- TRUE  # Set TRUE to enable combo analysis

# ==== HELPER FUNCTION ====
process_acn_group <- function(df, group_label, group_type) {
  df <- df %>%
    filter(
      !is.na(PURPL_ACN), PURPL_ACN >= 2, PURPL_ACN <= 6,
      !is.na(overallSurvival), !is.na(deceased)
    ) %>%
    mutate(
      overallSurvival = as.numeric(as.character(overallSurvival)),
      deceased = as.numeric(as.character(deceased))
    ) %>%
    filter(!is.na(overallSurvival), !is.na(deceased))
  #Check for variation in ACN
  
  if (nrow(df) < 2 || length(unique(df$deceased)) < 2) {
    return(data.frame(
      Group_Type = group_type,
      Group = group_label,
      n = nrow(df),
      PURPL_Cutoff = NA,
      MaxStat_statistic = NA,
      PURPL_p_value = NA
    ))
  }
  
  cutpoint <- tryCatch({
    surv_cutpoint(df, time = "overallSurvival", event = "deceased", variables = "PURPL_ACN", minprop = 0.1)
  }, error = function(e) NULL)
  
  if (is.null(cutpoint) || is.infinite(cutpoint$cutpoint[[1]]) || is.na(cutpoint$cutpoint[[1]])) {
    return(data.frame(
      Group_Type = group_type,
      Group = group_label,
      n = nrow(df),
      PURPL_Cutoff = NA,
      MaxStat_statistic = NA,
      PURPL_p_value = NA
    ))
  }
  
  cutoff_val <- cutpoint$cutpoint[[1]]
  maxstat_stat <- summary(cutpoint)$statistic[[1]]
  
  df <- df %>%
    mutate(PURPL_Group = ifelse(PURPL_ACN <= cutoff_val, "Low", "High")) %>%
    filter(!is.na(PURPL_Group))
  
  if (length(unique(df$PURPL_Group)) < 2) {
    return(data.frame(
      Group_Type = group_type,
      Group = group_label,
      n = nrow(df),
      PURPL_Cutoff = round(cutoff_val, 3),
      MaxStat_statistic = round(maxstat_stat, 4),
      PURPL_p_value = NA
    ))
  }
  
  pval <- tryCatch({
    surv_obj <- Surv(df$overallSurvival, df$deceased)
    fit <- survdiff(surv_obj ~ PURPL_Group, data = df)
    p <- pchisq(fit$chisq, df = 1, lower.tail = FALSE)
    signif(p, 4)
  }, error = function(e) NA)
  
  return(data.frame(
    Group_Type = group_type,
    Group = group_label,
    n = nrow(df),
    PURPL_Cutoff = round(cutoff_val, 3),
    MaxStat_statistic = round(maxstat_stat, 4),
    PURPL_p_value = pval
  ))
}

# ==== PART 1: DRIVER + MUTATION COMBINATIONS ====
acn_combo_results <- list()

if (combine_statuses) {
  combo_gene_pairs <- list(
    c("PIK3CA", "TP53"),
    c("PIK3CA", "AKT1")
  )
  
  for (pair in combo_gene_pairs) {
    gene1 <- pair[1]
    gene2 <- pair[2]
    
    col1 <- paste0(gene1, "_Mut_Status")
    col2 <- paste0(gene2, "_Mut_Status")
    
    if (!all(c(col1, col2) %in% colnames(TCGA_geneOfInterest))) next
    
    for (status1 in c(paste0(gene1, "_WT"), paste0(gene1, "_Mut"))) {
      for (status2 in c(paste0(gene2, "_WT"), paste0(gene2, "_Mut"))) {
        
        df_sub <- TCGA_geneOfInterest %>%
          filter(!!sym(col1) == status1, !!sym(col2) == status2)
        
        label <- paste(status1, "+", status2)
        
        result <- process_acn_group(df_sub, label, paste0(gene1, " + ", gene2, " Combo"))
        acn_combo_results[[length(acn_combo_results) + 1]] <- result
      }
    }
  }
}

# ==== PART 2: SINGLE MUTATION STATUS ====
acn_single_results <- list()

mutation_genes1 <- unique(c(mutation_genes, sub("_.*", "", driver_statuses)))
statuses_to_test <- unlist(lapply(mutation_genes1, function(g) c(paste0(g, "_WT"), paste0(g, "_Mut"))))

for (status in statuses_to_test) {
  gene <- sub("_.*", "", status)
  col <- paste0(gene, "_Mut_Status")
  if (!col %in% colnames(TCGA_geneOfInterest)) next
  
  df_sub <- TCGA_geneOfInterest %>% filter(!!sym(col) == status)
  result <- process_acn_group(df_sub, status, "Single Mutation Status")
  acn_single_results[[length(acn_single_results) + 1]] <- result
}

# ==== PART 3: GLOBAL ACN SURVIVAL ====
global_result <- process_acn_group(
  TCGA_geneOfInterest,
  group_label = "All Samples",
  group_type = "Overall Expression Effect"
)

# ==== PART 4: CUSTOM TRIPLE-MUTATION GROUPS ====
acn_triple_combo_results <- list()

status_cols <- c("PIK3CA_Mut_Status", "AKT1_Mut_Status", "MTOR_Mut_Status")

if (all(status_cols %in% colnames(TCGA_geneOfInterest))) {
  triple_groups <- list(
    list(
      label = "PIK3CA_WT + AKT1_WT + MTOR_WT",
      filters = c("PIK3CA_WT", "AKT1_WT", "MTOR_WT")
    ),
    list(
      label = "PIK3CA_WT + AKT1_WT + MTOR_Mut",
      filters = c("PIK3CA_WT", "AKT1_WT", "MTOR_Mut")
    )
  )
  
  for (group_def in triple_groups) {
    group_label <- group_def$label
    filter_vals <- group_def$filters
    
    df_sub <- TCGA_geneOfInterest %>%
      filter(
        PIK3CA_Mut_Status == filter_vals[1],
        AKT1_Mut_Status == filter_vals[2],
        MTOR_Mut_Status == filter_vals[3]
      )
    
    result <- process_acn_group(df_sub, group_label, "PIK3CA + AKT1 + MTOR Combo")
    if (!is.null(result)) acn_triple_combo_results[[length(acn_triple_combo_results) + 1]] <- result
  }
}

# ==== COMBINE RESULTS ====
acn_final <- bind_rows(acn_combo_results, acn_single_results, acn_triple_combo_results)

# ==== P-VALUE ADJUSTMENT ====
acn_final <- acn_final %>%
  mutate(Adjusted_p = p.adjust(PURPL_p_value, method = "BH"))

# ==== ADD GLOBAL RESULT ====
if (!is.null(global_result)) {
  acn_final <- bind_rows(acn_final, global_result)
}

# ==== OUTPUT ====
print(acn_final)
#write.csv(acn_final, "~/Data/TCGA/TCGA_pval_acn.csv", row.names = FALSE)
################################################################################





################################################################################
##Expression Survival using MaxStat adjusted p-value (BH) Truncation 2500 days##

# == Helper function to apply administrative censoring (truncation 2500 days) ==
cap_followup <- function(dat, time_col = "overallSurvival", event_col = "deceased", tau = 2500) {
  time  <- suppressWarnings(as.numeric(as.character(dat[[time_col]])))
  event <- suppressWarnings(as.numeric(as.character(dat[[event_col]])))
  ok <- !is.na(time) & !is.na(event)
  
  # Replace follow-up times exceeding tau with tau, and mark as censored
  dat[[time_col]][ok]  <- pmin(time[ok], tau)
  dat[[event_col]][ok & time[ok] > tau] <- 0
  dat
}


# == Data subset: Define the comparison group by uncommenting one or more mutation 
# filters below. This uses the corresponding adjusted p-value from `expr_tbl`.==

df_pik3ca <- TCGA_geneOfInterest %>%
  filter(
    #TP53_Mut_Status == "TP53_WT",
    PIK3CA_Mut_Status == "PIK3CA_WT",
    #AKT1_Mut_Status == "AKT1_WT",
    #ER_Status == "Negative",
    #HER2_Status == "Negative",
    !is.na(PURPL_Exp),
    !is.na(PURPL_ACN),
    PURPL_ACN >= 2,     # Lower bound of acceptable copy number
    PURPL_ACN <= 6,     # Upper bound of acceptable copy number
    !is.na(overallSurvival),
    !is.na(deceased)
  ) %>%
  mutate(
    overallSurvival = as.numeric(overallSurvival),
    deceased        = as.numeric(deceased)
  )

# Apply censoring
df_pik3ca <- cap_followup(df_pik3ca, "overallSurvival", "deceased", tau = 2500)


# ===== MaxStat cutoff for PURPL expression ======
cutpoint <- tryCatch(
  surv_cutpoint(df_pik3ca, time = "overallSurvival", event = "deceased", variables = "PURPL_Exp"),
  error = function(e) NULL
)

if (is.null(cutpoint)) {
  plot.new()
  title("MaxStat could not determine a cutoff for PURPL_Exp")
} else {
  # Create high/low groups using the MaxStat cutoff
  df_pik3ca <- df_pik3ca %>%
    mutate(PURPL_HighLow = ifelse(PURPL_Exp <= cutpoint$cutpoint[[1]], "Low", "High"))
  
  # Compute Kaplan–Meier fit and raw log-rank p-value
  surv_obj <- Surv(df_pik3ca$overallSurvival, df_pik3ca$deceased)
  fit <- survfit(surv_obj ~ PURPL_HighLow, data = df_pik3ca)
  raw_p <- as.numeric(surv_pvalue(fit, data = df_pik3ca)$pval)
  
  mut_cols <- c("PIK3CA_Mut_Status","TP53_Mut_Status","AKT1_Mut_Status","MTOR_Mut_Status")
  mut_cols <- mut_cols[mut_cols %in% colnames(df_pik3ca)]
  
  # Determine which genes have fixed values in the filtered data
  fixed <- lapply(mut_cols, function(cl) {
    vals <- unique(df_pik3ca[[cl]])
    vals <- vals[!is.na(vals)]
    if (length(vals) == 1) as.character(vals) else NA_character_
  })
  names(fixed) <- mut_cols
  fixed <- Filter(Negate(is.na), fixed)
  
  # Build candidate lookup keys for expr_tbl
  lookup_rows <- NULL
  if (length(fixed) >= 2) {
    # Two fixed genes → search for both orderings in expr_tbl
    genes <- sub("_Mut_Status$","", names(fixed))
    g1 <- genes[1]; g2 <- genes[2]
    s1 <- fixed[[1]]; s2 <- fixed[[2]]
    group_candidates <- c(paste0(s1, " + ", s2), paste0(s2, " + ", s1))
    type_candidates  <- c(paste0(g1, " + ", g2, " Combo"), paste0(g2, " + ", g1, " Combo"))
    lookup_rows <- subset(expr_tbl, Group %in% group_candidates & Group_Type %in% type_candidates)
    
  } else if (length(fixed) == 1) {
    # One fixed gene → single mutation status
    s1 <- unname(fixed[[1]])
    lookup_rows <- subset(expr_tbl, Group_Type == "Single Mutation Status" & Group == s1)
    
  } else {
    # No fixed genes → all samples (i.e. PURPL gene background)
    lookup_rows <- subset(expr_tbl, Group_Type == "Overall Expression Effect" & Group == "All Samples")
  }
  
  disp_p <- raw_p
  if (nrow(lookup_rows) >= 1) {
    p_bh  <- suppressWarnings(as.numeric(lookup_rows$Adjusted_p[1]))
    p_raw <- suppressWarnings(as.numeric(lookup_rows$PURPL_p_value[1]))
    if (!is.na(p_bh))      disp_p <- p_bh
    else if (!is.na(p_raw)) disp_p <- p_raw
  }
  pval_text <- paste0("p = ", signif(disp_p, 4))

  # ====== Plot: Kaplan–Meier with risk table (styled) =======
  g <- ggsurvplot(
    fit,
    data = df_pik3ca,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.title = NULL,
    conf.int = FALSE,
    legend = c(0.85, 0.35),
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
  
  # Remove redundant axis labels from the risk table
  g$table <- g$table +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  print(g)
}
################################################################################





################################################################################
#### ACN Survival using MaxStat adjusted p-value (BH) Truncation 8000 days #####

# === Fuction to ap follow-up at 8000 days" ===
cap_followup <- function(dat, time_col = "overallSurvival", event_col = "deceased", tau = 2500) {
  # "Administrative censoring: truncate times > tau and set event = 0 there"
  time  <- suppressWarnings(as.numeric(as.character(dat[[time_col]])))
  event <- suppressWarnings(as.numeric(as.character(dat[[event_col]])))
  ok <- !is.na(time) & !is.na(event)
  dat[[time_col]][ok]  <- pmin(time[ok], tau)
  dat[[event_col]][ok & time[ok] > tau] <- 0
  dat
}

# === "Filter & preprocess (uncomment what you need)" ===
df_km_base <- TCGA_geneOfInterest %>%
  filter(
    #PIK3CA_Mut_Status == "PIK3CA_WT",
    #TP53_Mut_Status   == "TP53_WT",
    #AKT1_Mut_Status   == "AKT1_WT",
    !is.na(PURPL_ACN),
    PURPL_ACN >= 2, PURPL_ACN <= 6,     # "Copy-number window"
    !is.na(overallSurvival), !is.na(deceased)
  ) %>%
  mutate(
    overallSurvival = as.numeric(overallSurvival),
    deceased        = as.numeric(deceased)
  ) %>%
  cap_followup("overallSurvival", "deceased", tau = 8000)

# === "MaxStat cutoff on PURPL_ACN" ===
cutpoint <- tryCatch(
  surv_cutpoint(df_km_base, time = "overallSurvival", event = "deceased",
                variables = "PURPL_ACN", minprop = 0.10),
  error = function(e) NULL
)

if (is.null(cutpoint) || is.na(cutpoint$cutpoint[[1]]) || is.infinite(cutpoint$cutpoint[[1]])) {
  plot.new(); title("MaxStat could not determine a cutoff for PURPL_ACN")
} else {
  optimal_cutoff <- cutpoint$cutpoint[[1]]
  
  # === "Create ACN strata (≤ cutoff vs > cutoff)" ===
  low_label  <- paste0("ACN \u2264 ", signif(optimal_cutoff, 4))
  high_label <- paste0("ACN > ",  signif(optimal_cutoff, 4))
  
  df_km <- df_km_base %>%
    mutate(
      PURPL_Group = ifelse(PURPL_ACN <= optimal_cutoff, low_label, high_label),
      PURPL_Group = factor(PURPL_Group, levels = c(high_label, low_label))  # "High first in legend"
    )
  
  # === "KM fit + raw log-rank p" ===
  surv_obj <- Surv(df_km$overallSurvival, df_km$deceased)
  fit <- survfit(surv_obj ~ PURPL_Group, data = df_km)
  lr  <- survdiff(surv_obj ~ PURPL_Group, data = df_km)
  pval_raw <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  
  # === "Figure out which BH-adjusted p to display (from acn_final)" ===
  # "Detect active mutation-status filters in the current df"
  mut_cols <- c("PIK3CA_Mut_Status","TP53_Mut_Status","AKT1_Mut_Status","MTOR_Mut_Status")
  mut_cols <- mut_cols[mut_cols %in% colnames(df_km)]
  
  fixed <- lapply(mut_cols, function(cl) {
    u <- unique(df_km[[cl]]); u <- u[!is.na(u)]
    if (length(u) == 1) as.character(u) else NA_character_
  })
  names(fixed) <- mut_cols
  fixed <- Filter(Negate(is.na), fixed)
  
  # "Build lookup against acn_final according to how many filters are fixed"
  lookup_rows <- NULL
  if (length(fixed) >= 2) {
    # "Two fixed genes → Combo row; accept either order for Group and Group_Type"
    genes <- sub("_Mut_Status$","", names(fixed))
    g1 <- genes[1]; g2 <- genes[2]
    s1 <- fixed[[1]]; s2 <- fixed[[2]]
    group_candidates <- c(paste0(s1, " + ", s2), paste0(s2, " + ", s1))
    type_candidates  <- c(paste0(g1, " + ", g2, " Combo"), paste0(g2, " + ", g1, " Combo"))
    lookup_rows <- subset(acn_final, Group %in% group_candidates & Group_Type %in% type_candidates)
    
  } else if (length(fixed) == 1) {
    # "One fixed gene → Single Mutation Status row"
    s1 <- unname(fixed[[1]])
    lookup_rows <- subset(acn_final, Group_Type == "Single Mutation Status" & Group == s1)
    
  } else {
    # "No fixed genes → Overall Expression Effect / All Samples row"
    lookup_rows <- subset(acn_final, Group_Type == "Overall Expression Effect" & Group == "All Samples")
  }
  
  # === "Prefer BH-adjusted p from acn_final; else raw p" ===
  disp_p <- pval_raw
  if (is.data.frame(lookup_rows) && nrow(lookup_rows) >= 1) {
    p_bh  <- suppressWarnings(as.numeric(lookup_rows$Adjusted_p[1]))
    p_raw <- suppressWarnings(as.numeric(lookup_rows$ACN_p_value[1] %||% NA_real_))
    if (!is.na(p_bh))      disp_p <- p_bh
    else if (!is.na(p_raw)) disp_p <- p_raw
  }
  pval_text <- paste0("p = ", signif(disp_p, 4))
  
  # === "Plot: KM with risk table (styled)" ===
  g <- ggsurvplot(
    fit,
    data = df_km,
    pval = pval_text,
    pval.size = 8,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.title = NULL,
    legend = c(0.85, 0.85),
    legend.title = "PURPL copy num.",
    legend.labs  = c(high_label, low_label),
    legend.title.font = 2,
    legend.text.font  = 2,
    xlab = "Time (days)",
    ylab = "Survival probability",
    font.x = c(30, "bold"),
    font.y = c(30, "bold"),
    font.tickslab = 30,
    font.legend = c(20, "bold"),
    tables.theme = theme_classic(base_size = 25),
    ggtheme = theme_classic() +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color = "black", fill = "transparent", size = 0.5)
      )
  )
  
  # === "Polish risk table: remove redundant axes" ===
  g$table <- g$table +
    labs(x = NULL, y = NULL, title = NULL) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  
  print(g)
}
################################################################################



