# Running guide to recurrent event discrimination on four models from SANAD and OPCRD datasets each

#=== setup ===#
library(Rcpp)
library(survival)
library(data.table)
library(ggplot2)
library(dplyr)
install.packages("tictoc")
library(tictoc)
library(MASS)

# Define the folder path
folder_path <- "M:/R/Data/Evaluating (calibration) paper/Predicted and observed datasets for both"
rds_files <- list.files(path = folder_path, pattern = "\\.rds$", full.names = TRUE)
rds_list <- lapply(rds_files, readRDS)
names(rds_list) <- tools::file_path_sans_ext(basename(rds_files))
list2env(rds_list, envir = .GlobalEnv)

# Cleaning these eight datasets
# Rename harmoniously
# nb_opcrd <- newdata_nb_opcrd
# zinb_opcrd <- newdata_zinb_opcrd
# ag_opcrd <- ag_pred_single_opcrd
# pwp_opcrd <- pwp_pred_single_opcrd
# nb_sanad <- newdata1c_nb_sanad
# zinb_sanad <- fullzdnbc_zinb_sanad
# ag_sanad <- newdataAGcfinal_ag_sanad
# pwp_sanad <- newdata4cfinal_pwp_sanad

# rm(newdata_nb_opcrd)
# rm(newdata_zinb_opcrd)
# rm(ag_pred_single_opcrd)
# rm(pwp_pred_single_opcrd)
# rm(newdata1c_nb_sanad)
# rm(fullzdnbc_zinb_sanad)
# rm(newdataAGcfinal_ag_sanad)
# rm(newdata4cfinal_pwp_sanad)

# List of all datasets
datasets <- list(
  nb_opcrd   = newdata_nb_opcrd,
  zinb_opcrd = newdata_zinb_opcrd,
  ag_opcrd   = ag_pred_single_opcrd,
  pwp_opcrd  = pwp_pred_single_opcrd,
  nb_sanad   = newdata1c_nb_sanad,
  zinb_sanad = fullzdnbc_zinb_sanad_v2,
  ag_sanad   = newdataAGcfinal_ag_sanad,
  pwp_sanad  = newdata4cfinal_pwp_sanad,
  sgf_sanad   = newdataSGFcfinal_SGF_sanad
)

# Clean each dataset
for (name in names(datasets)) {
  df <- datasets[[name]]

  # Coerce to data.frame early (before subsetting or renaming)
  df <- as.data.frame(df)
  
  # Rename ID variable
  if ("Patient_ID" %in% names(df)) {
    names(df)[names(df) == "Patient_ID"] <- "ID"
  } else if ("trialno" %in% names(df)) {
    names(df)[names(df) == "trialno"] <- "ID"
  }

  # Sanad-specific cleaning
  if (grepl("sanad", name)) {
    if ("totsez" %in% names(df)) names(df)[names(df) == "totsez"] <- "ObsCount"
    if ("TotSez" %in% names(df)) names(df)[names(df) == "TotSez"] <- "ObsCount"
    if ("SeizureCount" %in% names(df)) names(df)[names(df) == "SeizureCount"] <- "PredCount"
    
  } else if (grepl("opcrd", name)) {
    if (all(c("AttackCount", "ExacCount") %in% names(df))) {
      df$AttackCount <- NULL
    }
    if ("AttackCount" %in% names(df)) {
      names(df)[names(df) == "AttackCount"] <- "PredCount"
    } else if ("ExacCount" %in% names(df)) {
      names(df)[names(df) == "ExacCount"] <- "PredCount"
    }
    if ("csumexac" %in% names(df)) {
      names(df)[names(df) == "csumexac"] <- "ObsCount"
    }
  }

  # Ensure only the relevant variables are kept, and in the right format
  keep_vars <- c("ID", "ObsCount", "PredCount")
  present_vars <- intersect(keep_vars, names(df))
  df <- df[, present_vars, drop = FALSE]
  df <- as.data.frame(df)  # just to be extra safe

  # Reassign cleaned dataframe back to original name
  assign(name, df, envir = .GlobalEnv)
}


# List of final dataset names to keep
final_datasets <- c(
  "nb_opcrd", "zinb_opcrd", "ag_opcrd", "pwp_opcrd",
  "nb_sanad", "zinb_sanad", "ag_sanad", "pwp_sanad"
)

# Remove all objects not in final_datasets
rm(list = setdiff(ls(), final_datasets))




#====### Discrimination ###====#
setwd("M:/R")
sourceCpp("jackknifing_C_index_v7.cpp")

# Grouped jackknifing as a function

#' Random grouped jackknife (fast, vectors API; R-side set generation)
#'
#' @param ids   Character (or coercible) vector of IDs (length n)
#' @param obs   Numeric vector of observed counts (length n)
#' @param pred  Numeric vector of predicted counts (length n)
#' @param type  Metric: "C type1","C type2","C type3","kendall","goodman","C type4","somer"
#' @param d     Number of observations to drop per replicate, or a percent if d_is_percent=TRUE
#' @param n_reps Number of replicates (rows of drop sets)
#' @param d_is_percent Logical; if TRUE, interpret d as a percentage (1..99)
#' @param seed  Optional RNG seed for reproducibility
#' @param return_sets Logical; if TRUE, return list(reps=..., drop_sets=...)
#'
#' @return Numeric vector of replicate estimates (attributes preserved from C++), or a list if return_sets=TRUE
gjk_random_overlapping_fast_vec <- function(ids, obs, pred,
                                            type,
                                            d, n_reps,
                                            d_is_percent = FALSE,
                                            seed = NULL,
                                            return_sets = FALSE,
                                            drop_sets = NULL) {
  # --- Coerce & validate inputs ---
  ids  <- as.character(ids)
  obs  <- as.numeric(obs)
  pred <- as.numeric(pred)

  n <- length(ids)
  if (length(obs) != n || length(pred) != n) stop("ids, obs, and pred must have the same length.")
  if (n < 2L) stop("Need at least 2 observations.")
  if (any(!is.finite(obs)) || any(!is.finite(pred))) stop("obs/pred contain NA/NaN/Inf — please clean these.")

  # Validate metric early (C++ also validates)
  valid_types <- c("C type1","C type2","C type3","kendall","goodman","C type4","somer")
  if (!type %in% valid_types) {
    stop(sprintf("Invalid 'type'. Choose one of: %s", paste(valid_types, collapse = ", ")))
  }

  if (!is.null(seed)) set.seed(seed)
  if (n_reps < 2L) stop("n_reps must be at least 2.")

  # --- Either accept prebuilt sets or construct them here ---
  if (!is.null(drop_sets)) {
    # Basic checks on provided sets
    drop_sets <- as.matrix(drop_sets)
    if (!is.integer(drop_sets)) storage.mode(drop_sets) <- "integer"
    if (ncol(drop_sets) <= 0L || ncol(drop_sets) >= n) stop("ncol(drop_sets) must be between 1 and n-1.")
    if (any(drop_sets < 1L | drop_sets > n)) stop("drop_sets contains out-of-range indices.")
    if (nrow(drop_sets) < 2L) stop("drop_sets must have at least 2 rows (replicates).")
  } else {
    # --- Handle d as percent or absolute ---
    if (isTRUE(d_is_percent)) {
      if (d <= 0 || d >= 100) stop("When d_is_percent=TRUE, d must be in 1..99 (percent of n).")
      d <- max(1, min(n - 1L, round((d / 100) * n)))
    } else {
      if (d <= 0 || d >= n) stop("d must be between 1 and n-1.")
    }
    d <- as.integer(d)

    # --- Build n_reps × d matrix of 1-based indices (each row = one replicate) ---
    # Explicit replace = FALSE for clarity
    drop_sets <- t(replicate(n_reps, sample.int(n, d, replace = FALSE), simplify = TRUE))
    storage.mode(drop_sets) <- "integer"
  }

  # --- Call the fast C++ implementation that accepts precomputed sets ---
  reps <- gjk_random_overlapping_fast_with_sets(
    ids       = ids,
    obs       = obs,
    pred      = pred,
    type      = type,
    drop_sets = drop_sets
  )

  if (isTRUE(return_sets)) {
    return(list(reps = reps, drop_sets = drop_sets))
  }
  reps
}


## Make benchmarks single-threaded where possible
install.packages("RhpcBLASctl")
library(RhpcBLASctl)

if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}
if (requireNamespace("data.table", quietly = TRUE)) {
  data.table::setDTthreads(1)
}


## Core specs ---------------------------------------------------------------
get_cpu <- function() {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Windows") {
    x <- try(system("powershell -Command \"(Get-CimInstance Win32_Processor).Name\"", intern = TRUE), silent = TRUE)
    if (inherits(x, "try-error") || length(x) == 0) x <- system("wmic cpu get name", intern = TRUE)
    x <- x[nzchar(x)]
    return(trimws(tail(x, 1)))
  } else if (sys == "Darwin") {
    return(system("/usr/sbin/sysctl -n machdep.cpu.brand_string", intern = TRUE))
  } else { # Linux/Unix
    x <- try(system("lscpu | grep 'Model name' | sed 's/Model name:\\s*//'", intern = TRUE), silent = TRUE)
    if (!inherits(x, "try-error") && length(x)) return(trimws(x))
    # fallback
    return("Unknown CPU")
  }
}

get_ram_gb <- function() {
  sys <- Sys.info()[["sysname"]]
  bytes <- NA_real_
  if (sys == "Windows") {
    x <- suppressWarnings(system("wmic ComputerSystem get TotalPhysicalMemory", intern = TRUE))
    x <- as.numeric(gsub("\\D", "", x))
    bytes <- max(x, na.rm = TRUE)
  } else if (sys == "Darwin") {
    x <- suppressWarnings(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)))
    bytes <- x
  } else { # Linux/Unix
    x <- try(readLines("/proc/meminfo"), silent = TRUE)
    if (!inherits(x, "try-error")) {
      mt <- sub("MemTotal:\\s+([0-9]+) kB", "\\1", grep("^MemTotal", x, value = TRUE))
      bytes <- as.numeric(mt) * 1024
    }
  }
  round(bytes / (1024^3), 1) # GB
}

spec <- list(
  CPU              = get_cpu(),
  Cores_logical    = parallel::detectCores(logical = TRUE),
  Cores_physical   = tryCatch(parallel::detectCores(logical = FALSE), error = function(e) NA_integer_),
  RAM_GB           = get_ram_gb(),
  OS               = paste(Sys.info()[c("sysname","release")], collapse = " "),
  R                = R.version.string,
  R_platform       = R.version$platform,
  Rcpp_version     = as.character(utils::packageVersion("Rcpp"))
)

spec

# Minimal Wald CI from grouped-jackknife replicates (partition-style formula)
gjk_wald_ci_simple <- function(reps, level = 0.95) {
  # Point estimate: prefer full-sample theta if attached; else mean of reps
  theta_hat <- attr(reps, "theta_full")
  if (!is.numeric(theta_hat) || !is.finite(theta_hat)) theta_hat <- mean(reps, na.rm = TRUE)

  G <- length(reps)
  theta_bar <- mean(reps, na.rm = TRUE)

  # Jackknife variance (partition formula)
  var_hat <- (G - 1) / G * sum((reps - theta_bar)^2, na.rm = TRUE)
  se <- sqrt(var_hat)

  z <- qnorm(1 - (1 - level) / 2)
  c(lower = theta_hat - z * se,
    upper = theta_hat + z * se)
}

## START LOOP HERE ##

# OPCRD takes too long to run so run just on SANAD for now. 30/07/25 now both run in an acceptable timeframe!
# Define dataset names
dataset_names <- c(
  "nb_opcrd", "zinb_opcrd", "ag_opcrd", "pwp_opcrd",
  "nb_sanad", "zinb_sanad", "ag_sanad", "pwp_sanad"
)

# Create a storage list
results_list <- list()

# Define function for standard error
se <- function(x) sd(x) / sqrt(length(x))


## -------- Single-run Grouped JK (d = 10%, 200 reps) with partition-Wald CI --------

# Settings
gjk_d_pct   <- 10
gjk_n_reps  <- 200
gjk_seed    <- 123
gjk_types   <- c("kendall", "goodman", "somer", "C type4")

# Wald CI using grouped-JK replicate vector (partition variance; same as gjk_wald_ci_simple)
gjk_ci_wald_partition_from_reps <- function(reps, level = 0.95) {
  reps <- as.numeric(reps)
  reps <- reps[is.finite(reps)]
  G <- length(reps)
  if (G < 2L) stop("Need at least 2 grouped-JK replicates.")

  theta_hat <- attr(reps, "theta_full")
  if (!is.numeric(theta_hat) || !is.finite(theta_hat)) theta_hat <- mean(reps, na.rm = TRUE)

  theta_bar <- mean(reps, na.rm = TRUE)
  var_hat   <- (G - 1) / G * sum((reps - theta_bar)^2, na.rm = TRUE)
  se        <- sqrt(var_hat)

  z <- qnorm(1 - (1 - level)/2)
  c(lower = theta_hat - z * se,
    upper = theta_hat + z * se,
    se    = se,
    width = 2 * z * se)
}

gjk10_results <- list()

for (name in dataset_names) {
  if (!exists(name, inherits = FALSE)) {
    warning(sprintf("Object '%s' not found; skipping.", name))
    next
  }
  message("\n[Grouped-JK 10% | 200 reps] Dataset: ", name)
  df <- get(name)

  for (tp in gjk_types) {
    message("  - Type: ", tp)

    # One randomized grouped-JK run (overlapping sets), 200 replicates, d = 10%
    t0 <- proc.time()[["elapsed"]]
    reps <- gjk_random_overlapping_fast_vec(
      ids  = as.character(df$ID),
      obs  = as.numeric(df$ObsCount),
      pred = as.numeric(df$PredCount),
      type = tp,
      d = gjk_d_pct,
      n_reps = gjk_n_reps,
      d_is_percent = TRUE,
      seed = gjk_seed,
      return_sets = FALSE
    )
    t1 <- proc.time()[["elapsed"]]
    elapsed_gjk <- t1 - t0

    G <- length(reps)

    # SE-based (Wald) CI using partition variance (as in gjk_wald_ci_simple)
    ci_wald <- gjk_ci_wald_partition_from_reps(reps, level = 0.95)

    # Percentile CI (unchanged logic)
    qs <- as.numeric(quantile(reps, c(0.025, 0.975), names = FALSE))

    # Robust summaries (unchanged logic)
    med_reps <- median(reps)
    mad_reps <- mad(reps)

    # Center estimate: prefer attr theta_full if present
    est_center <- attr(reps, "theta_full")
    if (!is.numeric(est_center) || !is.finite(est_center)) {
      est_center <- mean(reps, na.rm = TRUE) # fallback
    }

    gjk10_results[[paste(name, tp, sep = "::")]] <- data.frame(
      dataset       = name,
      type          = tp,
      d_label       = "10%",
      G             = G,
      est           = est_center,
      # SE-based (partition) CI
      ci_wald_lower = ci_wald["lower"],
      ci_wald_upper = ci_wald["upper"],
      se_wald       = ci_wald["se"],
      width_wald    = ci_wald["width"],
      # Percentile CI (unchanged)
      ci_pct_lower  = qs[1],
      ci_pct_upper  = qs[2],
      width_pct     = qs[2] - qs[1],
      # Robust replicate summaries (unchanged)
      median_reps   = med_reps,
      mad_reps      = mad_reps,
      # Timing
      elapsed_gjk   = elapsed_gjk,
      stringsAsFactors = FALSE
    )
  }
}

gjk10_df <- do.call(rbind, gjk10_results)

# Optional compact summary
gjk10_summary <- gjk10_df %>%
  dplyr::group_by(dataset, type, d_label) %>%
  dplyr::summarise(
    G                = unique(G),
    est              = unique(est),
    width_wald       = unique(width_wald),
    width_pct        = unique(width_pct),
    median_reps      = unique(median_reps),
    mad_reps         = unique(mad_reps),
    elapsed_gjk      = unique(elapsed_gjk),
    .groups = "drop"
  )

options(scipen = 999, digits = 15)
View(gjk10_df)
View(gjk10_summary)

write.csv(gjk10_df,      "GroupedJK_10pct_200reps_runlevel.csv", row.names = FALSE)
write.csv(gjk10_summary, "GroupedJK_10pct_200reps_summary.csv",  row.names = FALSE)

