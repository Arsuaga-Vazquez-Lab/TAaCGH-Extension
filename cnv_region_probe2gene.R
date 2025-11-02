# ===========================================================
# End-to-end: probe stats -> significance -> boxplots -> genes -> patient profiles
# ===========================================================
# ---- Packages ----
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("biomaRt")
}
library(dplyr)
library(biomaRt)

# ---- Modify as needed ----
phenColumn <- "Luminal_A"  # e.g. "PIK3CA_mut", "Luminal_A", "Luminal_B", "Basal", "HER2", etc.
chr        <- "5"          # chromosome without "chr"
arm        <- "p"          # "p" or "q"
startBP    <- 19608109     # region start bp
endBP      <- 42861888     # region end bp
cghStart   <- 6            # first sample column in "_data_full.txt"

# ---- File pathway, please modify as necessary ----
dataRoot    <- path.expand("~/Desktop/P1BCB215/Data/Horlings")
resultsRoot <- path.expand("~/Desktop/P1BCB215/Results/Horlings")
dataPath    <- file.path(dataRoot, "horlings_data_full.txt")
phenPath    <- file.path(dataRoot, "horlings_phen.txt")

# ---- Load data ----
data     <- read.table(dataPath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
phenData <- read.table(phenPath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!(phenColumn %in% colnames(phenData))) stop(paste("Column", phenColumn, "not found in phenData"))


# ================================
# Region slice
# ================================
regionData <- data %>%
  dplyr::filter(Chrom == chr, Arm == arm, bp >= startBP, bp <= endBP) %>%
  dplyr::arrange(bp)
if (nrow(regionData) == 0) stop("No probes found in the requested region.")
copyMat <- regionData[, cghStart:ncol(regionData), drop = FALSE]
rownames(copyMat) <- regionData$Clone

# ================================
# Phenotype alignment
# ================================
sample_ids_phen  <- as.character(phenData[, 1])    # assuming col 1 is sample ID
sample_ids_cgh   <- colnames(copyMat)
sample_ids_phen_prefixed <- paste0("X", sample_ids_phen)

matched <- match(sample_ids_cgh, sample_ids_phen_prefixed)
phenVec <- rep(NA, length(sample_ids_cgh))
phenVec[!is.na(matched)] <- phenData[matched[!is.na(matched)], phenColumn]

if (sum(!is.na(phenVec)) < length(phenVec) / 3) {
  matched2 <- match(sample_ids_cgh, sample_ids_phen)
  phenVec2 <- rep(NA, length(sample_ids_cgh))
  phenVec2[!is.na(matched2)] <- phenData[matched2[!is.na(matched2)], phenColumn]
  if (sum(!is.na(phenVec2)) > sum(!is.na(phenVec))) phenVec <- phenVec2
}
if (sum(!is.na(phenVec)) == 0) stop("No matching phenotype data for CGH columns.")

make_group_labels <- function(phen_col) {
  if (grepl("mut", phen_col, ignore.case = TRUE)) {
    list(g1 = phen_col, g0 = sub("(?i)mut", "WT", phen_col, perl = TRUE))
  } else {
    list(g1 = paste0(phen_col, "_pos"), g0 = paste0(phen_col, "_neg"))
  }
}
glabs <- make_group_labels(phenColumn)

# ================================
# Permutation test
# ================================
permutation_test <- function(group1, group0, n_perm = 100000, stat_func = function(x, y) mean(x) - mean(y)) {
  obs_stat <- stat_func(group1, group0)
  combined <- c(group1, group0); n1 <- length(group1)
  perm_stats <- replicate(n_perm, {
    perm <- sample(combined)
    stat_func(perm[1:n1], perm[(n1 + 1):length(perm)])
  })
  p_value <- mean(abs(perm_stats) >= abs(obs_stat))
  c(p_value = p_value, stat = obs_stat)
}

# ================================
# Probe-wise stats
# ================================
results <- data.frame()

for (i in seq_len(nrow(copyMat))) {
  probeVals <- suppressWarnings(as.numeric(copyMat[i, ]))
  probeID   <- regionData$Clone[i]
  bp        <- regionData$bp[i]
  
  valid <- !is.na(phenVec) & !is.na(probeVals)
  if (sum(valid) < 4) next
  
  groupNum <- ifelse(phenVec[valid] == 1, 1, ifelse(phenVec[valid] == 0, 0, NA))
  keep <- !is.na(groupNum)
  if (sum(keep) < 4) next
  
  df <- data.frame(copyNumber = probeVals[valid][keep],
                   groupNum   = factor(groupNum[keep], levels = c(0, 1)))
  
  if (length(unique(df$groupNum)) == 2 && all(table(df$groupNum) >= 2)) {
    wilcox_res <- suppressWarnings(wilcox.test(copyNumber ~ groupNum, data = df, exact = FALSE))
    ttest_res  <- suppressWarnings(t.test(copyNumber ~ groupNum, data = df))
    perm_res   <- permutation_test(df$copyNumber[df$groupNum == 1],
                                   df$copyNumber[df$groupNum == 0])
    
    mean_g1 <- mean(df$copyNumber[df$groupNum == 1])
    mean_g0 <- mean(df$copyNumber[df$groupNum == 0])
    delta   <- mean_g1 - mean_g0
    direction <- ifelse(delta > 0, "Amplification",
                        ifelse(delta < 0, "Deletion", "NoChange"))
    
    results <- rbind(results, data.frame(
      Probe       = probeID,
      Chr         = chr,
      Arm         = arm,
      bp          = bp,
      Mean_Group0 = mean_g0,
      Mean_Group1 = mean_g1,
      Delta_G1_G0 = delta,
      Direction   = direction,
      Wilcox_p    = wilcox_res$p.value,
      Ttest_p     = ttest_res$p.value,
      Ttest_stat  = as.numeric(ttest_res$statistic),
      Perm_p      = perm_res["p_value"],
      Perm_stat   = perm_res["stat"],
      stringsAsFactors = FALSE
    ))
  }
}
if (nrow(results) == 0) stop("No analyzable probes (insufficient samples per group).")

# ================================
# Multiple testing (Bonferroni and BH)
# ================================
results <- results %>%
  dplyr::mutate(
    Bonferroni_Wilcox = p.adjust(Wilcox_p, method = "bonferroni"),
    Bonferroni_Ttest  = p.adjust(Ttest_p,  method = "bonferroni"),
    Bonferroni_Perm   = p.adjust(Perm_p,   method = "bonferroni"),
    BH_Wilcox         = p.adjust(Wilcox_p, method = "BH"),
    BH_Ttest          = p.adjust(Ttest_p,  method = "BH"),
    BH_Perm           = p.adjust(Perm_p,   method = "BH")
  )

# ================================
# Output folders
# ================================
mb_start <- formatC(startBP / 1e6, format = "f", digits = 2)
mb_end   <- formatC(endBP   / 1e6, format = "f", digits = 2)
region_folder_name <- paste0("chr", chr, arm, "_", mb_start, "MB-", mb_end, "MB")

baseOutDir <- file.path(resultsRoot, phenColumn, region_folder_name)
if (!dir.exists(baseOutDir)) dir.create(baseOutDir, recursive = TRUE)

# Master stats CSV for the region
master_stats_path <- file.path(baseOutDir, paste0("chr", chr, arm, "_", startBP, "_", endBP, "_", phenColumn, "_stats.csv"))
write.csv(results, file = master_stats_path, row.names = FALSE)

# ================================
# Patient profiles (phenotype==1)
# ================================
phen_pos_idx <- which(phenVec == 1 & !is.na(phenVec))
if (length(phen_pos_idx) > 0) {
  Mbp  <- regionData$bp / 1e6
  subMat <- suppressWarnings(apply(copyMat[, phen_pos_idx, drop = FALSE], 2, as.numeric))
  if (is.vector(subMat)) subMat <- matrix(subMat, ncol = 1)
  y_min <- suppressWarnings(min(subMat, na.rm = TRUE)); y_max <- suppressWarnings(max(subMat, na.rm = TRUE))
  if (!is.finite(y_min) || !is.finite(y_max)) { y_min <- -1; y_max <- 1 }
  pad   <- 0.05 * (y_max - y_min)
  ylim  <- c(y_min - pad, y_max + pad)
  
  page_title <- sprintf("%s – chr%s%s, %s–%s Mb",
                        phenColumn, chr, arm,
                        formatC(startBP/1e6, format = "f", digits = 2),
                        formatC(endBP/1e6,   format = "f", digits = 2))
  
  prof_pdf <- file.path(baseOutDir,
                        paste0(phenColumn, "_patient_profiles_", "chr", chr, arm, "_",
                               formatC(startBP/1e6, format="f", digits=2), "MB-",
                               formatC(endBP/1e6,   format="f", digits=2), "MB.pdf"))
  
  pdf(prof_pdf, width = 9, height = 6)
  plots_per_page <- 6
  par(mfrow = c(3, 2),
      mar  = c(4, 4, 2.2, 1) + 0.1,   # inner margins
      oma  = c(0, 0, 2.5, 0))         # outer margin for page title
  
  for (j in seq_along(phen_pos_idx)) {
    col_idx <- phen_pos_idx[j]
    pat_id  <- colnames(copyMat)[col_idx]
    pat_vals <- suppressWarnings(as.numeric(copyMat[, col_idx]))
    
    plot(Mbp, pat_vals,
         type = "o", pch = 20, col = "gray20",
         main = paste("Patient", pat_id),
         xlab = "Position (Mb)",
         ylab = "log2 copy-number ratio",
         ylim = ylim)
    abline(h = 0, lty = 2, col = "black")
    
    # If we've filled a page or reached the last patient, write the page title
    if ((j %% plots_per_page == 0) || (j == length(phen_pos_idx))) {
      mtext(page_title, outer = TRUE, cex = 1.1, line = 0.5, font = 2)
      # If more plots remain, start a new page with the same layout
      if (j < length(phen_pos_idx)) {
        par(mfrow = c(3, 2),
            mar  = c(4, 4, 2.2, 1) + 0.1,
            oma  = c(0, 0, 2.5, 0))
      }
    }
  }
  dev.off()
  message("Saved patient profiles PDF: ", prof_pdf)
} else {
  message("No patients with phenotype '", phenColumn, "' == 1 for profiles. Skipping PDF.")
}

# ================================
# Significant probes
# ================================
sig_any <- results %>%
  dplyr::mutate(
    Sig_Wilcox_Bonf = Bonferroni_Wilcox < 0.05,
    Sig_Ttest_Bonf  = Bonferroni_Ttest  < 0.05,
    Sig_Perm_Bonf   = Bonferroni_Perm   < 0.05,
    Sig_Wilcox_BH   = BH_Wilcox         < 0.05,
    Sig_Ttest_BH    = BH_Ttest          < 0.05,
    Sig_Perm_BH     = BH_Perm           < 0.05,
    Bonferroni = Sig_Wilcox_Bonf | Sig_Ttest_Bonf | Sig_Perm_Bonf,
    BH         = Sig_Wilcox_BH   | Sig_Ttest_BH   | Sig_Perm_BH
  ) %>%
  dplyr::arrange(dplyr::desc(Bonferroni), dplyr::desc(BH),
                 Perm_p, Wilcox_p, Ttest_p)

sig_probes_tbl <- sig_any %>%
  dplyr::filter(Bonferroni | BH) %>%
  dplyr::select(Probe, bp, Direction, Bonferroni, BH)

sig_probes_path <- file.path(baseOutDir, "significant_probes.csv")
write.csv(sig_probes_tbl, file = sig_probes_path, row.names = FALSE)

sig_probes <- sig_probes_tbl$Probe
if (length(sig_probes) == 0) {
  message("No significant probes by any corrected test. No per-probe folders will be created.")
}

# ================================
# biomaRt gene fetcher
# ================================
get_genes_from_ensembl <- function(chr_label_with_chr, start, end, version = 109) {
  ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", version = version)
  res <- getBM(
    attributes = c("hgnc_symbol","ensembl_gene_id","gene_biotype",
                   "chromosome_name","start_position","end_position","strand"),
    filters = c("chromosome_name","start","end"),
    values = list(sub("^chr", "", chr_label_with_chr), start, end),
    mart = ensembl
  )
  if (nrow(res) == 0) return(res)
  res[order(res$start_position), , drop = FALSE]
}

# Adjacent interval helper
adjacent_interval_for_probe <- function(probe_id) {
  i <- which(regionData$Clone == probe_id)
  if (length(i) != 1) return(NULL)
  left_bp  <- if (i > 1) regionData$bp[i - 1] else regionData$bp[i]
  right_bp <- if (i < nrow(regionData)) regionData$bp[i + 1] else regionData$bp[i]
  c(min(left_bp, right_bp), max(left_bp, right_bp))
}

# ================================
# Per-probe outputs
# - PDF boxplot (line1: probe/bp/direction; line2: p-values in smaller font)
# - boxplot_stats.csv (merged sections)
# - genes_between_adjacent_probes.csv with appended summaries + overlapping genes
# ================================
for (p in sig_probes) {
  subdir <- file.path(baseOutDir, p)
  if (!dir.exists(subdir)) dir.create(subdir, recursive = TRUE)
  
  file.copy(master_stats_path, file.path(subdir, "all_probes_significance.csv"), overwrite = TRUE)
  
  i <- which(regionData$Clone == p)
  if (length(i) != 1) next
  probe_bp <- regionData$bp[i]
  probeVals <- suppressWarnings(as.numeric(copyMat[i, ]))
  valid     <- !is.na(phenVec) & !is.na(probeVals)
  gnum <- ifelse(phenVec[valid] == 1, 1, ifelse(phenVec[valid] == 0, 0, NA))
  keep <- !is.na(gnum)
  
  dfBox <- data.frame(
    copyNumber = probeVals[valid][keep],
    groupNum   = factor(gnum[keep], levels = c(0, 1))
  )
  dfBox$group <- factor(ifelse(dfBox$groupNum == 1, glabs$g1, glabs$g0),
                        levels = c(glabs$g0, glabs$g1))
  
  # Pull p-values & direction(amplification/deletion) from results
  rrow   <- results[results$Probe == p, ][1, ]
  dir_txt <- rrow$Direction
  p_w <- signif(rrow$Wilcox_p, 3)
  p_t <- signif(rrow$Ttest_p,  3)
  
  # ---- PDF boxplot ----
  y_min <- suppressWarnings(min(dfBox$copyNumber, na.rm = TRUE))
  y_max <- suppressWarnings(max(dfBox$copyNumber, na.rm = TRUE))
  y_rng <- if (is.finite(y_min) && is.finite(y_max)) (y_max - y_min) else 1
  ylim  <- c(y_min, y_max + 0.15 * y_rng)
  
  pdf(file.path(subdir, paste0(p, "_", phenColumn, "_boxplot.pdf")), width = 6.5, height = 5.8)
  par(mar = c(4.1, 4.1, 4.2, 2.1))  # extra top margin
  
  main_line1 <- sprintf("Probe %s, bp=%s, %s", p, format(probe_bp, big.mark=","), dir_txt)
  boxplot(
    copyNumber ~ group,
    data    = dfBox,
    main    = main_line1,
    ylab    = "log2 copy-number ratio",
    xlab    = "",
    outline = TRUE,
    col     = c("lightgreen", "skyblue"),
    ylim    = ylim,
    cex.main = 0.95
  )
  abline(h = 0, lty = 2)
  # second line with p-values
  mtext(sprintf("p=%s (Wilcoxon), p=%s (T-test)", p_w, p_t),
        side = 3, line = 0.6, cex = 0.8)
  
  dev.off()
  
  # --- Combined boxplot stats into a single CSV ---
  dfClean <- dfBox[!is.na(dfBox$group) & !is.na(dfBox$copyNumber), ]
  dfClean$group <- droplevels(dfClean$group)
  
  summary_stats <- dfClean %>%
    dplyr::group_by(group) %>%
    dplyr::summarize(
      N      = dplyr::n(),
      Median = median(copyNumber, na.rm = TRUE),
      Q1     = quantile(copyNumber, 0.25, na.rm = TRUE),
      Q3     = quantile(copyNumber, 0.75, na.rm = TRUE),
      Mean   = mean(copyNumber, na.rm = TRUE),
      SD     = sd(copyNumber, na.rm = TRUE),
      .groups = "drop"
    )
  
  stats_out <- file.path(subdir, "boxplot_stats.csv")
  writeLines("summary_by_group", con = stats_out)
  suppressWarnings(write.table(summary_stats, file = stats_out, sep = ",",
                               row.names = FALSE, col.names = TRUE, quote = TRUE, append = TRUE))
  
  two_group_ok <- (nlevels(dfClean$group) == 2) && all(table(dfClean$group) >= 2)
  wilcox_p <- NA_real_
  ttest_p  <- NA_real_
  if (two_group_ok) {
    wilcox_p <- suppressWarnings(wilcox.test(copyNumber ~ group, data = dfClean, exact = FALSE)$p.value)
    ttest_p  <- suppressWarnings(t.test(copyNumber ~ group, data = dfClean)$p.value)
  }
  cat("\n\ntwo_sample_tests\n", file = stats_out, append = TRUE)
  suppressWarnings(write.table(data.frame(Wilcoxon_p = wilcox_p, Ttest_p = ttest_p),
                               file = stats_out, sep = ",", row.names = FALSE,
                               col.names = TRUE, quote = TRUE, append = TRUE))
  
  do_vs0 <- function(vals) {
    data.frame(
      wilcox_p = suppressWarnings(wilcox.test(vals, mu = 0)$p.value),
      ttest_p  = suppressWarnings(t.test(vals,  mu = 0)$p.value)
    )
  }
  vs0_df <- data.frame(Group = character(), wilcox_p = numeric(), ttest_p = numeric())
  wt_level <- levels(dfClean$group)[grepl("WT", levels(dfClean$group), ignore.case = TRUE)][1]
  if (!is.na(wt_level)) {
    tmp <- do_vs0(dfClean$copyNumber[dfClean$group == wt_level]); tmp$Group <- wt_level
    vs0_df <- tmp[, c("Group","wilcox_p","ttest_p")]
  } else {
    for (lvl in levels(dfClean$group)) {
      tmp <- do_vs0(dfClean$copyNumber[dfClean$group == lvl]); tmp$Group <- lvl
      vs0_df <- rbind(vs0_df, tmp[, c("Group","wilcox_p","ttest_p")])
    }
  }
  cat("\n\none_sample_vs0\n", file = stats_out, append = TRUE)
  suppressWarnings(write.table(vs0_df, file = stats_out, sep = ",",
                               row.names = FALSE, col.names = TRUE, quote = TRUE, append = TRUE))
  
  # ---- Genes between adjacent probes + summaries + overlapping genes ----
  iv <- adjacent_interval_for_probe(p)
  if (!is.null(iv)) {
    genes_df <- get_genes_from_ensembl(chr_label_with_chr = paste0("chr", chr),
                                       start = iv[1], end = iv[2], version = 109)
    if (nrow(genes_df) > 0) {
      genes_df$overlaps_probe_bp <- (genes_df$start_position <= probe_bp) & (probe_bp <= genes_df$end_position)
    } else {
      genes_df$overlaps_probe_bp <- logical(0)
    }
    
    genes_out <- file.path(subdir, "genes_between_adjacent_probes.csv")
    write.csv(genes_df, file = genes_out, row.names = FALSE)
    
    if (nrow(genes_df) > 0) {
      biotype_counts <- as.data.frame(sort(table(genes_df$gene_biotype), decreasing = TRUE))
      colnames(biotype_counts) <- c("gene_biotype", "count")
    } else {
      biotype_counts <- data.frame(gene_biotype = character(0), count = integer(0))
    }
    cat("\n\nbiotype_summary\n", file = genes_out, append = TRUE)
    suppressWarnings(write.table(biotype_counts, file = genes_out, sep = ",",
                                 row.names = FALSE, col.names = TRUE, quote = TRUE, append = TRUE))
    
    overlapping <- genes_df[genes_df$overlaps_probe_bp, , drop = FALSE]
    cat("\n\ngenes_overlapping_probe_bp\n", file = genes_out, append = TRUE)
    if (nrow(overlapping) > 0) {
      suppressWarnings(write.table(overlapping[, c("hgnc_symbol","ensembl_gene_id","gene_biotype",
                                                   "chromosome_name","start_position","end_position","strand")],
                                   file = genes_out, sep = ",", row.names = FALSE,
                                   col.names = TRUE, quote = TRUE, append = TRUE))
    } else {
      suppressWarnings(write.table(data.frame(hgnc_symbol=character(0),
                                              ensembl_gene_id=character(0),
                                              gene_biotype=character(0),
                                              chromosome_name=character(0),
                                              start_position=integer(0),
                                              end_position=integer(0),
                                              strand=integer(0)),
                                   file = genes_out, sep = ",", row.names = FALSE,
                                   col.names = TRUE, quote = TRUE, append = TRUE))
    }
  }
}

cat("\n--- OUTPUT ---\n")
cat("Region folder: ", baseOutDir, "\n")
cat("Master stats:  ", master_stats_path, "\n")
cat("Significant:   ", sig_probes_path, "\n")
cat("Per-probe folders include: boxplot PDF (2-line title), COMBINED boxplot_stats.csv, genes CSV with appended summaries and overlapping genes.\n")

