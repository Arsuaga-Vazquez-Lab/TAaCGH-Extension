# ============================
# METABRIC PURPL CNA → KM Plot
# ============================

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# ---- File paths (edit if needed) ----
clinical_path <- "/Users/salvador/Downloads/brca_metabric_clinical_data.tsv"
cna_path      <- "/Users/salvador/Downloads/METABRICCNA_chr5.csv"

# ---- PURPL locus (hg19) ----
#purpl_start <- 27472292L
#purpl_end   <- 27496400L

# ---- PURPL locus (hg38) ----
#purpl_start <- 27472246L
#purpl_end   <- 27496404L

purpl_start <- 27217714L
purpl_end   <- 27497871L

# ---- Helpers ----
norm_names <- function(nm) {
  nm |>
    tolower() |>
    str_replace_all("[^a-z0-9]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")
}

pick_col <- function(df, options) {
  opts <- norm_names(options)
  names(df) <- norm_names(names(df))
  hit <- opts[opts %in% names(df)]
  if (!length(hit)) stop("None of these columns found: ", paste(options, collapse = ", "))
  hit[1]
}

# =========================
# 1) Clinical → survival
# =========================
clinical_raw <- read_tsv(clinical_path, show_col_types = FALSE)
names(clinical_raw) <- norm_names(names(clinical_raw))

id_col    <- pick_col(clinical_raw, c("patient_id","metabric_id","sample_id","id"))
os_m_col  <- pick_col(clinical_raw, c("overall_survival_months","os_months"))
os_s_col  <- pick_col(clinical_raw, c("overall_survival_status","os_status"))
pam50_col <- pick_col(clinical_raw, c("pam50_claudin_low_subtype","pam50","pam50_subtype"))

clinical_surv <- clinical_raw %>%
  dplyr::mutate(
    METABRIC_ID   = .data[[id_col]],
    pam50_subtype = .data[[pam50_col]],
    OS_months     = suppressWarnings(as.numeric(.data[[os_m_col]])),
    overallSurvival = OS_months * 30.44,  # months → days
    deceased = dplyr::case_when(
      is.na(.data[[os_s_col]]) ~ NA_integer_,
      stringr::str_detect(.data[[os_s_col]], stringr::regex("^\\s*1\\b|deceased|dead", ignore_case = TRUE)) ~ 1L,
      stringr::str_detect(.data[[os_s_col]], stringr::regex("^\\s*0\\b|alive|living",  ignore_case = TRUE)) ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  dplyr::filter(tolower(pam50_subtype) == "luma") %>%
  dplyr::transmute(METABRIC_ID, overallSurvival, deceased)

# Quick sanity
cat("LumA clinical rows:", nrow(clinical_surv), "\n")

# =========================
# 2) CNA → PURPL status
# =========================
cna_raw <- read_csv(cna_path, show_col_types = FALSE)
names(cna_raw) <- norm_names(names(cna_raw))

# Resolve columns (robust to case/header variants)
cna_id    <- pick_col(cna_raw, c("metabric_id","patient_id","sample_id","id"))
chrom_col <- pick_col(cna_raw, c("chrom","chromosome"))
start_col <- pick_col(cna_raw, c("loc_start","start"))
end_col   <- pick_col(cna_raw, c("loc_end","end"))
call_col  <- pick_col(cna_raw, c("call2","call","segment_call"))

cna_data <- cna_raw %>%
  dplyr::transmute(
    METABRIC_ID = .data[[cna_id]],
    chrom       = .data[[chrom_col]],
    loc_start   = suppressWarnings(as.numeric(.data[[start_col]])),
    loc_end     = suppressWarnings(as.numeric(.data[[end_col]])),
    call2       = toupper(as.character(.data[[call_col]]))
  )

# Counts
cat("Unique METABRIC_IDs in CNA file:", dplyr::n_distinct(cna_data$METABRIC_ID), "\n")
cat("Unique METABRIC_IDs in clinical_surv:", dplyr::n_distinct(clinical_surv$METABRIC_ID), "\n")
cat("Shared IDs:", length(intersect(unique(cna_data$METABRIC_ID), unique(clinical_surv$METABRIC_ID))), "\n")

# PURPL status (GAIN/AMP vs NEUT) based on overlap with locus window
purpl_status <- cna_data %>%
  dplyr::filter(chrom %in% c(5, "5")) %>%
  dplyr::filter(loc_end >= purpl_start, loc_start <= purpl_end) %>%
  dplyr::filter(call2 %in% c("AMP", "GAIN", "NEUT")) %>%
  dplyr::mutate(PURPL_status = dplyr::case_when(
    call2 %in% c("AMP", "GAIN") ~ "GAIN_OR_AMP",
    call2 == "NEUT" ~ "NEUT"
  )) %>%
  dplyr::group_by(METABRIC_ID) %>%
  # Prefer GAIN/AMP over NEUT when both exist
  dplyr::slice_max(
    order_by = factor(PURPL_status, levels = c("NEUT", "GAIN_OR_AMP")),
    n = 1, with_ties = FALSE
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(METABRIC_ID, PURPL_status)

# =========================
# 3) Merge & KM
# =========================
merged_data <- clinical_surv %>%
  dplyr::inner_join(purpl_status, by = "METABRIC_ID") %>%
  dplyr::mutate(
    PURPL_status = factor(PURPL_status, levels = c("GAIN_OR_AMP", "NEUT"))
  )

cat("PURPL status table:\n")
print(table(merged_data$PURPL_status))

# KM fit + log-rank p
surv_obj <- survival::Surv(time = merged_data$overallSurvival, event = merged_data$deceased)
fit <- survival::survfit(surv_obj ~ PURPL_status, data = merged_data)

lr  <- survival::survdiff(surv_obj ~ PURPL_status, data = merged_data)
p_raw <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
formatted_pval <- sprintf("p = %.4f", p_raw)

# =========================
# 4) Plot (your preferred style; no xlim/max-time)
# =========================
g <- survminer::ggsurvplot(
  fit,
  data = merged_data,
  
  pval = formatted_pval,
  pval.size = 8,
  
  conf.int = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  
  # keep colored strata tags in risk table
  risk.table.y.text = TRUE,
  risk.table.y.text.col = TRUE,
  risk.table.title = NULL,  # remove "Number at risk"
  
  legend = c(0.85, 0.85),
  legend.title = "PURPL copy num.",
  legend.labs  = c("GAIN", "NEUT"),
  legend.title.font = 2,
  legend.text.font  = 2,
  
  xlab = "Time (days)",
  ylab = "Survival probability",
  
  # big fonts
  font.x = c(30, "bold"),
  font.y = c(30, "bold"),
  font.tickslab = 30,
  font.legend = c(20, "bold"),
  
  tables.theme = ggplot2::theme_classic(base_size = 25),
  ggtheme = ggplot2::theme_classic() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(color = "black", fill = "transparent", size = 0.5)
    )
)

# Remove axis titles from risk table; keep colored tags
g$table <- g$table +
  ggplot2::labs(x = NULL, y = NULL, title = NULL) +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank()
  )

print(g)




























































suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

# ---- File paths (edit if needed) ----
clinical_path <- "/Users/salvador/Downloads/brca_metabric_clinical_data.tsv"
cna_path      <- "/Users/salvador/Downloads/METABRICCNA_chr5.csv"

# ---- PURPL locus (choose the genome build that matches your CNA) ----
# hg19:
# purpl_start <- 27472292L
# purpl_end   <- 27496400L
# hg38:
purpl_start <- 27217714L
purpl_end   <- 27497871L

# =========================
# Helpers
# =========================
norm_names <- function(nm) {
  nm |>
    tolower() |>
    stringr::str_replace_all("[^a-z0-9]+", "_") |>
    stringr::str_replace_all("_+", "_") |>
    stringr::str_replace_all("^_|_$", "")
}

pick_col <- function(df, options) {
  opts <- norm_names(options)
  names(df) <- norm_names(names(df))
  hit <- opts[opts %in% names(df)]
  if (!length(hit)) stop("None of these columns found: ", paste(options, collapse = ", "))
  hit[1]
}

# Administrative censoring at a chosen tau (days)
cap_followup <- function(dat, time_col = "overallSurvival", event_col = "deceased", tau) {
  time  <- suppressWarnings(as.numeric(as.character(dat[[time_col]])))
  event <- suppressWarnings(as.numeric(as.character(dat[[event_col]])))
  ok <- !is.na(time) & !is.na(event)
  time_cap  <- time
  event_cap <- event
  time_cap[ok] <- pmin(time[ok], tau)
  event_cap[ok & time[ok] > tau] <- 0
  dat[[time_col]]  <- time_cap
  dat[[event_col]] <- event_cap
  dat
}

logrank_p <- function(dat, group_col = "PURPL_status") {
  fit <- survival::survdiff(Surv(dat$overallSurvival, dat$deceased) ~ dat[[group_col]])
  stats::pchisq(fit$chisq, df = 1, lower.tail = FALSE)
}

summarise_groups <- function(d) {
  d %>% dplyr::group_by(PURPL_status) %>%
    dplyr::summarise(n = dplyr::n(), events = sum(deceased == 1, na.rm = TRUE), .groups = "drop")
}

plot_km <- function(dat, pval, tau) {
  surv_obj <- Surv(dat$overallSurvival, dat$deceased)
  fit <- survfit(surv_obj ~ PURPL_status, data = dat)
  lbl <- sprintf("p = %s", ifelse(is.na(pval), "NA", sprintf("%.4f", pval)))
  survminer::ggsurvplot(
    fit, data = dat,
    pval = lbl, pval.size = 8,
    conf.int = FALSE,
    risk.table = TRUE, risk.table.height = 0.25,
    risk.table.y.text = TRUE, risk.table.y.text.col = TRUE,
    risk.table.title = NULL,
    legend = c(0.85, 0.85),
    legend.title = "PURPL copy num.",
    legend.labs  = c("NEUT","GAIN"),
    legend.title.font = 2, legend.text.font = 2,
    xlab = "Time (days)", ylab = "Survival probability",
    xlim = c(0, tau),
    ggtheme = ggplot2::theme_classic() + ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.background = ggplot2::element_rect(color = "black", fill = "transparent", size = 0.5)
    ),
    tables.theme = ggplot2::theme_classic(base_size = 25)
  )
}

# =========================
# 1) Clinical → survival (LumA)
# =========================
clinical_raw <- readr::read_tsv(clinical_path, show_col_types = FALSE)
names(clinical_raw) <- norm_names(names(clinical_raw))

id_col    <- pick_col(clinical_raw, c("patient_id","metabric_id","sample_id","id"))
os_m_col  <- pick_col(clinical_raw, c("overall_survival_months","os_months"))
os_s_col  <- pick_col(clinical_raw, c("overall_survival_status","os_status"))
pam50_col <- pick_col(clinical_raw, c("pam50_claudin_low_subtype","pam50","pam50_subtype"))

clinical_surv <- clinical_raw %>%
  dplyr::mutate(
    METABRIC_ID     = .data[[id_col]],
    pam50_subtype   = .data[[pam50_col]],
    OS_months       = suppressWarnings(as.numeric(.data[[os_m_col]])),
    overallSurvival = OS_months * 30.44,  # months → days
    deceased        = dplyr::case_when(
      is.na(.data[[os_s_col]]) ~ NA_integer_,
      stringr::str_detect(.data[[os_s_col]], stringr::regex("^\\s*1\\b|deceased|dead", ignore_case = TRUE)) ~ 1L,
      stringr::str_detect(.data[[os_s_col]], stringr::regex("^\\s*0\\b|alive|living",  ignore_case = TRUE)) ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  dplyr::filter(tolower(pam50_subtype) == "luma") %>%
  dplyr::transmute(METABRIC_ID, overallSurvival, deceased)

cat("LumA clinical rows:", nrow(clinical_surv), "\n")

# =========================
# 2) CNA → PURPL status (GAIN/AMP vs NEUT)
# =========================
cna_raw <- readr::read_csv(cna_path, show_col_types = FALSE)
names(cna_raw) <- norm_names(names(cna_raw))

cna_id    <- pick_col(cna_raw, c("metabric_id","patient_id","sample_id","id"))
chrom_col <- pick_col(cna_raw, c("chrom","chromosome"))
start_col <- pick_col(cna_raw, c("loc_start","start","segment_start"))
end_col   <- pick_col(cna_raw, c("loc_end","end","segment_end"))
call_col  <- pick_col(cna_raw, c("call2","call","segment_call"))

cna_data <- cna_raw %>%
  dplyr::transmute(
    METABRIC_ID = .data[[cna_id]],
    chrom       = .data[[chrom_col]],
    loc_start   = suppressWarnings(as.numeric(.data[[start_col]])),
    loc_end     = suppressWarnings(as.numeric(.data[[end_col]])),
    call2       = toupper(as.character(.data[[call_col]]))
  )

cat("Unique METABRIC_IDs in CNA file:", dplyr::n_distinct(cna_data$METABRIC_ID), "\n")
cat("Unique METABRIC_IDs in clinical_surv:", dplyr::n_distinct(clinical_surv$METABRIC_ID), "\n")
cat("Shared IDs:", length(intersect(unique(cna_data$METABRIC_ID), unique(clinical_surv$METABRIC_ID))), "\n")

purpl_status <- cna_data %>%
  dplyr::filter(chrom %in% c(5, "5"),
                loc_end >= purpl_start, loc_start <= purpl_end) %>%
  dplyr::filter(call2 %in% c("AMP","GAIN","NEUT")) %>%
  dplyr::mutate(PURPL_status = dplyr::case_when(
    call2 %in% c("AMP","GAIN") ~ "GAIN",
    call2 == "NEUT"            ~ "NEUT"
  )) %>%
  dplyr::group_by(METABRIC_ID) %>%
  # Prefer GAIN over NEUT if both appear
  dplyr::slice_max(order_by = factor(PURPL_status, levels = c("NEUT","GAIN")),
                   n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(METABRIC_ID, PURPL_status)

# =========================
# 3) Merge & prepare analysis frame
# =========================
dat_all <- clinical_surv %>%
  dplyr::inner_join(purpl_status, by = "METABRIC_ID") %>%
  dplyr::filter(!is.na(overallSurvival), !is.na(deceased)) %>%
  dplyr::mutate(
    overallSurvival = suppressWarnings(as.numeric(overallSurvival)),
    deceased        = suppressWarnings(as.numeric(deceased)),
    PURPL_status    = factor(PURPL_status, levels = c("NEUT","GAIN"))
  ) %>%
  dplyr::filter(!is.na(overallSurvival), !is.na(deceased))

cat("PURPL status table (LumA):\n"); print(table(dat_all$PURPL_status))
if (nrow(dat_all) < 3 || length(unique(dat_all$deceased)) < 2) {
  stop("Not enough usable LumA rows after merge (need ≥3 and at least two event states).")
}

# =========================
# 4) Run analyses at τ = 2500 and τ = 8000 days
# =========================
taus <- c(2500, 8000)

results <- lapply(taus, function(tau) {
  d_tau <- cap_followup(dat_all, "overallSurvival", "deceased", tau = tau)
  p_tau <- if (sum(d_tau$deceased == 1, na.rm = TRUE) == 0 ||
               length(unique(d_tau$deceased)) < 2) NA_real_ else logrank_p(d_tau, "PURPL_status")
  s_tau <- summarise_groups(d_tau)
  list(tau = tau, dat = d_tau, p = p_tau, sum = s_tau)
})

# Build a tidy results table
row_for_tau <- function(res) {
  data.frame(
    tau_days            = res$tau,
    n_NEUT              = with(res$sum, n[PURPL_status=="NEUT"]),
    events_NEUT         = with(res$sum, events[PURPL_status=="NEUT"]),
    n_GAIN              = with(res$sum, n[PURPL_status=="GAIN"]),
    events_GAIN         = with(res$sum, events[PURPL_status=="GAIN"]),
    p_logrank           = signif(res$p, 4),
    check.names = FALSE
  )
}
res_table <- do.call(rbind, lapply(results, row_for_tau))

cat("\n=== METABRIC LumA — PURPL GAIN vs NEUT (administrative censoring) ===\n")
print(res_table, row.names = FALSE)

# =========================
# 5) Plots for each τ
# =========================
g_list <- mapply(function(res, tau) plot_km(res$dat, res$p, tau),
                 results, taus, SIMPLIFY = FALSE)

# Print plots (first τ=2500, then τ=8000)
for (g in g_list) print(g)

# =========================
# 5) Plots for each τ using your exact style
# =========================
make_plot_in_your_style <- function(dat_in, p_val, tau) {
  # Reorder levels so legend.labs = c("GAIN","NEUT") aligns correctly
  dat_plot <- dat_in
  dat_plot$PURPL_status <- factor(dat_plot$PURPL_status, levels = c("GAIN","NEUT"))
  
  # KM fit on the truncated (administratively censored) data for this tau
  surv_obj <- survival::Surv(dat_plot$overallSurvival, dat_plot$deceased)
  fit <- survival::survfit(surv_obj ~ PURPL_status, data = dat_plot)
  
  # Your p-value label format
  formatted_pval <- sprintf("p = %s", ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)))
  
  g <- survminer::ggsurvplot(
    fit,
    data = dat_plot,
    
    pval = formatted_pval,
    pval.size = 8,
    
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    
    # keep colored strata tags in risk table
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.title = NULL,  # remove "Number at risk"
    
    legend = c(0.85, 0.85),
    legend.title = "PURPL copy num.",
    legend.labs  = c("GAIN", "NEUT"),
    legend.title.font = 2,
    legend.text.font  = 2,
    
    xlab = "Time (days)",
    ylab = "Survival probability",
    
    # big fonts
    font.x = c(30, "bold"),
    font.y = c(30, "bold"),
    font.tickslab = 30,
    font.legend = c(20, "bold"),
    
    tables.theme = ggplot2::theme_classic(base_size = 25),
    ggtheme = ggplot2::theme_classic() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(color = "black", fill = "transparent", size = 0.5)
      )
  )
  
  # Remove axis titles from risk table; keep colored tags
  g$table <- g$table +
    ggplot2::labs(x = NULL, y = NULL, title = NULL) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  
  print(g)
}

# Render one plot per tau using the truncated datasets you already created
for (res in results) {
  make_plot_in_your_style(dat_in = res$dat, p_val = res$p, tau = res$tau)
}







