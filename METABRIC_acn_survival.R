################################################################################
#======== File paths (download files before running) ========
clinical_path <- "~/Data/METABRIC/brca_metabric_clinical_data.tsv"
cna_path      <- "~/Data/METABRIC/METABRICCNADiscovery.txt"   
################################################################################
#======== Please load the libraries before every use ======== 
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(readr)
  library(survival)
  library(survminer)
  library(ggplot2)
})
################################################################################
# Survival p values and adjusted p values (BH) using Exp. Truncation 2500 days #

#======== PURPL locus (pick matching build) ========
# If your METABRIC CNA is GRCh37/hg19, use:
# purpl_start <- 27472292L; purpl_end <- 27496400L
# Below is the range from hg38
purpl_start <- 27217714L
purpl_end   <- 27497871L

#======== Helper functions ========
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

cap_followup <- function(dat, time_col = "overallSurvival", event_col = "deceased", tau = 8000) {
  time  <- suppressWarnings(as.numeric(dat[[time_col]]))
  event <- suppressWarnings(as.numeric(dat[[event_col]]))
  ok <- !is.na(time) & !is.na(event)
  time_cap  <- pmin(time, tau)
  event_cap <- ifelse(time > tau, 0, event)
  dat[[time_col]]  <- time_cap
  dat[[event_col]] <- event_cap
  dat
}

logrank_p <- function(dat, group_col = "PURPL_status") {
  s <- survival::survdiff(Surv(dat$overallSurvival, dat$deceased) ~ dat[[group_col]])
  stats::pchisq(s$chisq, df = 1, lower.tail = FALSE)
}

#======== KM plotting function ========
plot_metabric_purpl_km <- function(dat_in, p_val, tau = 8000) {
  dat_in$PURPL_status <- factor(dat_in$PURPL_status, levels = c("GAIN","NEUT"))
  surv_obj <- survival::Surv(dat_in$overallSurvival, dat_in$deceased)
  fit <- survival::survfit(surv_obj ~ PURPL_status, data = dat_in)
  formatted_pval <- sprintf("p = %s", ifelse(is.na(p_val), "NA", sprintf("%.4f", p_val)))
  
  g <- survminer::ggsurvplot(
    fit,
    data = dat_in,
    pval = formatted_pval,
    pval.size = 8,
    conf.int = FALSE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text = TRUE,
    risk.table.y.text.col = TRUE,
    risk.table.title = NULL,
    legend = c(0.85, 0.85),
    legend.title = "PURPL copy num.",
    legend.labs  = c("GAIN", "NEUT"),
    legend.title.font = 2,
    legend.text.font  = 2,
    xlab = "Time (days)",
    ylab = "Survival probability",
    xlim = c(0, tau),
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
  
  g$table <- g$table +
    ggplot2::labs(x = NULL, y = NULL, title = NULL) +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )
  
  print(g)
}

#======== Load clinical data and derive LumA survival ========
clinical_raw <- read_tsv(clinical_path, show_col_types = FALSE)
names(clinical_raw) <- norm_names(names(clinical_raw))

id_col    <- pick_col(clinical_raw, c("patient_id","metabric_id","sample_id","id"))
os_m_col  <- pick_col(clinical_raw, c("overall_survival_months","os_months"))
os_s_col  <- pick_col(clinical_raw, c("overall_survival_status","os_status"))
pam50_col <- pick_col(clinical_raw, c("pam50_claudin_low_subtype","pam50","pam50_subtype"))

clinical_surv <- clinical_raw %>%
  mutate(
    METABRIC_ID     = .data[[id_col]],
    pam50_subtype   = .data[[pam50_col]],
    OS_months       = suppressWarnings(as.numeric(.data[[os_m_col]])),
    overallSurvival = OS_months * 30.44,
    deceased        = case_when(
      is.na(.data[[os_s_col]]) ~ NA_integer_,
      str_detect(.data[[os_s_col]], regex("^\\s*1\\b|deceased|dead", ignore_case = TRUE)) ~ 1L,
      str_detect(.data[[os_s_col]], regex("^\\s*0\\b|alive|living",  ignore_case = TRUE)) ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  filter(tolower(pam50_subtype) == "luma") %>%
  transmute(METABRIC_ID, overallSurvival, deceased)

#======== Load CNA file and identify PURPL region ========
cna_raw <- read_tsv(cna_path, show_col_types = FALSE, progress = FALSE)
names(cna_raw) <- norm_names(names(cna_raw))

cna_id    <- pick_col(cna_raw, c("metabric_id","sample_id","sample","patient_id","id"))
chrom_col <- pick_col(cna_raw, c("chrom","chromosome"))
start_col <- pick_col(cna_raw, c("loc_start","start","segment_start","start_position"))
end_col   <- pick_col(cna_raw, c("loc_end","end","segment_end","end_position"))

# Prefer call column if available; otherwise use segment mean
call_col <- intersect(norm_names(c("call2","call","segment_call","cna_call")), names(cna_raw))[1]

cna_data <- cna_raw %>%
  transmute(
    METABRIC_ID = .data[[cna_id]],
    chrom       = .data[[chrom_col]],
    loc_start   = suppressWarnings(as.numeric(.data[[start_col]])),
    loc_end     = suppressWarnings(as.numeric(.data[[end_col]])),
    call_raw    = if (!is.na(call_col)) toupper(as.character(.data[[call_col]])) else NA_character_
  )

cna_chr5 <- cna_data %>%
  filter(chrom %in% c("5", 5),
         loc_end >= purpl_start,
         loc_start <= purpl_end)

#======== Define PURPL GAIN vs NEUT ========
if (!is.na(call_col)) {
  purpl_status <- cna_chr5 %>%
    filter(call_raw %in% c("AMP","GAIN","NEUT")) %>%
    mutate(PURPL_status = case_when(
      call_raw %in% c("AMP","GAIN") ~ "GAIN",
      call_raw == "NEUT"            ~ "NEUT"
    )) %>%
    group_by(METABRIC_ID) %>%
    slice_max(order_by = factor(PURPL_status, levels = c("NEUT","GAIN")),
              n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(METABRIC_ID, PURPL_status)
} else {
  seg_col <- pick_col(cna_raw, c("seg_mean","segment_mean","log2ratio","value"))
  cna_chr5_seg <- cna_raw
  names(cna_chr5_seg) <- norm_names(names(cna_chr5_seg))
  cna_chr5_seg <- cna_chr5_seg %>%
    transmute(
      METABRIC_ID = .data[[cna_id]],
      chrom       = .data[[chrom_col]],
      loc_start   = suppressWarnings(as.numeric(.data[[start_col]])),
      loc_end     = suppressWarnings(as.numeric(.data[[end_col]])),
      seg_mean    = suppressWarnings(as.numeric(.data[[seg_col]]))
    ) %>%
    filter(chrom %in% c("5",5),
           loc_end >= purpl_start,
           loc_start <= purpl_end)
  
  purpl_status <- cna_chr5_seg %>%
    mutate(PURPL_status = ifelse(seg_mean >= 0.3, "GAIN", "NEUT")) %>%
    group_by(METABRIC_ID) %>%
    slice_max(order_by = factor(PURPL_status, levels = c("NEUT","GAIN")),
              n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(METABRIC_ID, PURPL_status)
}

#======== Merge datasets, censor at 8000 days, and plot ========
tau <- 8000

dat_all <- clinical_surv %>%
  inner_join(purpl_status, by = "METABRIC_ID") %>%
  filter(!is.na(overallSurvival), !is.na(deceased)) %>%
  mutate(
    overallSurvival = as.numeric(overallSurvival),
    deceased        = as.numeric(deceased),
    PURPL_status    = factor(PURPL_status, levels = c("GAIN","NEUT"))
  )

dat_tau <- cap_followup(dat_all, "overallSurvival", "deceased", tau = tau)

p_tau <- if (length(unique(dat_tau$deceased)) >= 2) logrank_p(dat_tau, "PURPL_status") else NA_real_

plot_metabric_purpl_km(dat_tau, p_tau, tau = tau)
################################################################################




