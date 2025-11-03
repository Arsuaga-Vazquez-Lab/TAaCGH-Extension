################################################################################
# SCAN-B (GSE96058) — preprocessing (PURPL in LumA)
# - Loads curated metadata + expression from the GSE96058 R package
# - Cleans/standardizes metadata and merges expression (incl. PURPL/LINC01021)
# - Optionally merges mutation calls from a user-provided TSV (e.g., Mutation Explorer export)
# - Produces a LumA Kaplan–Meier with MaxStat cutpoint on PURPL expression
# - Outputs data file into ~/Data/SCAN-B
################################################################################

################################################################################
#======== Install and load libraries ===========================================
pkgs <- c("GSE96058","dplyr","stringr","tidyr","readr","survival","survminer")
to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if ("GSE96058" %in% to_install) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("12379Monty/GSE96058")
  to_install <- setdiff(to_install, "GSE96058")
}
if (length(to_install)) install.packages(to_install)

suppressPackageStartupMessages({
  library(GSE96058)
  library(dplyr); library(stringr); library(tidyr); library(readr)
  library(survival); library(survminer)
})
################################################################################

################################################################################
#======== Set pathways ==================================================================
outdir <- path.expand("~/Data/SCAN-B")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
################################################################################

################################################################################
#======== Load metadata and expression =========================================
meta <- get("GSE96058_sampDesc")

obj_names  <- ls("package:GSE96058")
expr_names <- grep("^GSE96058_geneExpression_sub[1-5]$", obj_names, value = TRUE)
stopifnot(length(expr_names) > 0)
expr_list <- lapply(expr_names, function(nm) get(nm, asNamespace("GSE96058")))
expr_all  <- do.call(cbind, expr_list)  # genes x samples (GEO colnames)
cat("Expression matrix:", nrow(expr_all), "genes x", ncol(expr_all), "samples\n")
################################################################################

################################################################################
#======== Clean / standardize metadata =========================================
meta_clean <- meta %>%
  rename_with(~ .x |>
                tolower() |>
                str_trim() |>
                str_replace_all("[\\s/]+","_") |>
                str_replace_all("[:;()]+","") |>
                str_replace_all("\\.+","_")) %>%
  mutate(across(everything(), as.character))
################################################################################

################################################################################
#======== Build meta_std =======================================================
meta_std <- meta_clean %>%
  mutate(
    isrepl_log = tolower(isrepl) %in% c("true","t","1","yes"),
    SAMPLE_ID  = coalesce(geoacc, scan_b_external_id, title, sampno),
    pam50_subtype = pam50_subtype,
    ER_status     = er_status,
    PR_status     = pgr_status,
    HER2_status   = her2_status,
    OS_days_raw   = ovrallsurvdays,
    OS_event_raw  = ovrallsurvevent,
    overallSurvival = suppressWarnings(as.numeric(OS_days_raw)),
    overallSurvival = ifelse(
      is.na(overallSurvival),
      suppressWarnings(as.numeric(stringr::str_extract(OS_days_raw, "\\d+\\.?\\d*"))),
      overallSurvival
    ),
    deceased = dplyr::case_when(
      is.na(OS_event_raw) ~ NA_integer_,
      tolower(OS_event_raw) %in% c("1","deceased","dead","event","yes") ~ 1L,
      tolower(OS_event_raw) %in% c("0","alive","censored","no") ~ 0L,
      suppressWarnings(!is.na(as.numeric(OS_event_raw))) ~ as.integer(as.numeric(OS_event_raw) > 0),
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(!isrepl_log | is.na(isrepl_log)) %>%
  dplyr::select(
    SAMPLE_ID, overallSurvival, deceased,
    pam50_subtype, ER_status, PR_status, HER2_status,
    everything()
  )
################################################################################

################################################################################
#======== Align samples between metadata and expression ========================
id_overlap <- intersect(colnames(expr_all), meta_std$title)
if (length(id_overlap) == 0)
  stop("No overlap between expression colnames and metadata 'title'")
expr_all <- expr_all[, id_overlap, drop = FALSE]
meta_std <- meta_std %>% filter(title %in% id_overlap)
################################################################################

################################################################################
#======== Genes of interest -> *_Exp columns ===================================
genes_of_interest <- c("LINC01021","TP53","PIK3CA","MTOR","ULK1","AKT1","NFKB1","NFKBIA","CDKN1A")
expr_rownames_up <- toupper(rownames(expr_all))
matched_genes <- rownames(expr_all)[expr_rownames_up %in% toupper(genes_of_interest)]
cat("Found genes in matrix:", paste(matched_genes, collapse = ", "), "\n")

expr_subset <- expr_all[matched_genes, , drop = FALSE]
GENE_Exp <- as.data.frame(t(expr_subset))
GENE_Exp$title <- rownames(GENE_Exp)

display_names <- matched_genes
display_names[grepl("^LINC01021$", display_names, ignore.case = TRUE)] <- "PURPL"
display_names <- paste0(display_names, "_Exp")
colnames(GENE_Exp) <- c(display_names, "title")

meta_std <- meta_std %>% left_join(GENE_Exp, by = "title")

write.csv(meta_std, file.path(outdir, "SCAN-B_GSE96058.csv"), row.names = FALSE)
cat("Wrote:", file.path(outdir, "SCAN-B_GSE96058.csv"), "\n")
################################################################################





################################################################################
#======== Mutation data merge ==================================================
# We use the downloaded mutation file from the SCAN-B Mutation Explorer:
# https://oncogenomics.bmc.lu.se/MutationExplorer/
# and relabeled as "SCAN-B_Mutation.tsv" to be included in /Data/SCAN-B folder
################################################################################
scanb_file <- file.path(outdir, "SCAN-B_GSE96058.csv")
meta_std <- read.csv(scanb_file, stringsAsFactors = FALSE, check.names = FALSE)

meta_std <- meta_std %>%
  mutate(SCANB_ID = ifelse(!is.na(scan_b_external_id),
                           stringr::str_extract(scan_b_external_id, "S\\d{6}"),
                           NA))

mutation_file <- "~/Data/SCAN-B/SCAN-B_Mutation.tsv"
mut_raw <- read_tsv(mutation_file, comment = "#", show_col_types = FALSE)
mutations <- mut_raw %>% mutate(SCANB_ID = stringr::str_extract(SAMPLE, "S\\d{6}"))

genes_of_interest <- c("TP53","PIK3CA","MTOR","ULK1","AKT1","NFKB1","NFKBIA","CDKN1A")
genes_present <- intersect(genes_of_interest, unique(mutations$gene.symbol))

if (length(genes_present) > 0) {
  mut_matrix <- mutations %>%
    filter(gene.symbol %in% genes_present) %>%
    mutate(Mutated = TRUE) %>%
    select(SCANB_ID, gene.symbol, Mutated) %>%
    distinct() %>%
    pivot_wider(names_from = gene.symbol, values_from = Mutated, values_fill = FALSE) %>%
    mutate(across(all_of(genes_present),
                  ~ ifelse(.x, paste0(cur_column(), "_Mut"), paste0(cur_column(), "_WT")))) %>%
    rename_with(~ paste0(.x, "_Mut_Status"), all_of(genes_present))
  
  meta_std <- meta_std %>% left_join(mut_matrix, by = "SCANB_ID")
  
  for (g in genes_present) {
    colname <- paste0(g, "_Mut_Status")
    if (colname %in% names(meta_std))
      meta_std[[colname]] <- ifelse(is.na(meta_std[[colname]]),
                                    paste0(g, "_WT"),
                                    meta_std[[colname]])
  }
}

write.csv(meta_std, scanb_file, row.names = FALSE)
cat("Updated file with mutation status:", scanb_file, "\n")
################################################################################












