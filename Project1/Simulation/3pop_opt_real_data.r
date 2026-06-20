rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

#============================================================
# REAL-DATA-SAFE PI ESTIMATION SCRIPT
#============================================================
# Goal:
#   1. Estimate pi from posterior soft counts.
#   2. Apply one fixed calibrated correction rule.
#   3. NEVER use true_pi to choose the final estimate.
#   4. Use true_pi only after estimation for simulation comparison.
#
# Important:
#   - Pi_Final_Optimized_NoTrue is the practical final estimate.
#   - The calibrated lambdas must be chosen from previous simulation summaries
#     across all scenarios, then fixed before applying this script to real data.
#   - For real data, set use_true_pi_for_comparison <- FALSE.
#============================================================

#============================================================
# INPUTS
#============================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Missing argument: SLURM_ARRAY_TASK_ID")
}

jset <- as.numeric(args[1])

file.rjobs <- "/lustre10/scratch/yatah3/yatah3/simulation/project1/step3/"

inputs <- read.table(
  paste0(file.rjobs, "input2.txt"),
  header = TRUE,
  as.is = TRUE
)

scenarioIndex <- inputs[jset, 1]

dirBase <- "/lustre10/scratch/yatah3/yatah3/simulation/SimuGenotype/"

scenario_dirs <- list.dirs(
  dirBase,
  full.names = TRUE,
  recursive = FALSE
)

scenario_dirs <- scenario_dirs[
  grepl("sim_hsq", basename(scenario_dirs))
]

scenario_dirs <- sort(scenario_dirs)

sc_dir <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)

density_dir <- paste0(sc_dir, "/new_densityU_3pop/")

summary_file <- list.files(
  path = sc_dir,
  pattern = "^MaxLogSum_sim(10|[1-9])\\.csv$",
  full.names = TRUE
)

if (length(summary_file) == 0) {
  stop("No summary files found in: ", sc_dir)
}

#============================================================
# SETTINGS
#============================================================

rho_vec <- c(seq(0, 0.9, 0.1), 0.95)

penalty <- "RealmixLOG"
warmStart <- 1
Zscale <- 1
singleStart <- 1

pi_names <- c(
  "pi_0",
  "pi_1",
  "pi_2",
  "pi_3",
  "pi_12",
  "pi_13",
  "pi_23",
  "pi_123"
)

#============================================================
# TRUE PI: COMPARISON ONLY
#============================================================
# TRUE VALUES ARE NOT USED TO COMPUTE OR SELECT THE FINAL ESTIMATE.
# They are only used after Pi_Final_Optimized_NoTrue is already computed.
#
# For real data, set this to FALSE.
# For simulation reports, keep TRUE.

use_true_pi_for_comparison <- TRUE

if (use_true_pi_for_comparison) {

pi_table <- matrix(
  c(
    # Scenario 1: new three-trait / one-pop result
    0.901457, 0.003889, 0.003713, 0.003670,
    0.011032, 0.011187, 0.010867, 0.054185,

    # Scenario 2: old result row 1
    0.702615263, 0.009838396, 0.009655417, 0.009595751,
    0.037893735, 0.038119143, 0.037961357, 0.154320939,

    # Scenario 3: new three-trait / one-pop result
    0.902320, 0.004186, 0.004036, 0.004163,
    0.009523, 0.009365, 0.009381, 0.057026,

    # Scenario 4: new three-trait / one-pop result
    0.906529, 0.004399, 0.004108, 0.004101,
    0.005061, 0.005037, 0.005100, 0.065665,

    # Scenario 5: old result row 2
    0.703071383, 0.010010767, 0.009992204, 0.010055848,
    0.036878073, 0.037002710, 0.036786584, 0.156202432,

    # Scenario 6: old result row 3
    0.703870918, 0.011052947, 0.010847428, 0.010807650,
    0.034466207, 0.034227540, 0.034276599, 0.160450710,

    # Scenario 7: old result row 4
    0.707051823, 0.012457106, 0.012190594, 0.012271476,
    0.028051361, 0.028549911, 0.028404058, 0.171023670,

    # Scenario 8: old result row 5
    0.712032013, 0.013234101, 0.013133330, 0.013345479,
    0.021545026, 0.021421714, 0.021510552, 0.183777785
  ),
  nrow = 8,
  byrow = TRUE
)

  colnames(pi_table) <- pi_names

  true_pi <- pi_table[scenarioIndex, ]
  names(true_pi) <- pi_names

} else {

  true_pi <- rep(NA_real_, length(pi_names))
  names(true_pi) <- pi_names
}

#============================================================
# CALIBRATED PARAMETERS FOR FINAL REAL-DATA-SAFE ESTIMATE
#============================================================
# These are the only tuning parameters used to compute the final estimate.
# They must be selected from aggregate simulation calibration across all scenarios.
#
# Based on your final simulation/excel summaries, corrected_pi123_v2 was the
# useful structure. The example values below should be replaced if your
# aggregate 5-scenario summary gives a different robust median/mode.
#
# Suggested from current results:
#   lambda_single_calibrated = 0.65
#   lambda_zero_calibrated   = 0.05
#
# Important:
#   These values are fixed before running this script on real data.
#   The script does not select lambdas using true_pi.

lambda_single_calibrated <- 0.50
lambda_zero_calibrated   <- 0.00

lambda_single_directional_calibrated <- 5
lambda_shared_directional_calibrated <- 0.025

#============================================================
# DIAGNOSTIC PARAMETER GRIDS
#============================================================
# These are saved for sensitivity analysis only.
# They are not used to select the final estimate.

correction_grid <- expand.grid(
  lambda_single = seq(0, 0.95, by = 0.05),
  lambda_zero   = c(0, 0.005, 0.01, 0.02, 0.03, 0.05),
  stringsAsFactors = FALSE
)


directional_lambda_grid <- expand.grid(
  lambda_single = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 0.80, 1, 1.5, 2, 5),
  lambda_shared = c(0, 0.01, 0.025, 0.05, 0.075, 0.10, 0.25, 0.50, 1, 2, 5),
  stringsAsFactors = FALSE
)

gamma_grid <- c(0.5, 1, 2, 3, 5, 10)

#============================================================
# HELPER FUNCTIONS
#============================================================

process_one_summary <- function(summary_file) {

  finfo <- file.info(summary_file)

  if (is.na(finfo$size) || finfo$size == 0) {
    message("Skipping empty file: ", summary_file)
    return(NULL)
  }

  best_res <- read.csv(summary_file, header = TRUE)

  if (nrow(best_res) == 0) {
    message("No rows in file: ", summary_file)
    return(NULL)
  }

  list(
    sim = best_res$sim[1],
    selected_roww = best_res$Largest_Index[1],
    selected_rho = best_res$Corresponding_Rho[1],
    tauuse = best_res$Corresponding_Tauuse[1],
    Largest_Sum_Of_Logs = best_res$Largest_Sum_Of_Logs[1]
  )
}

project_simplex <- function(x, eps = 1e-12) {

  x[!is.finite(x)] <- 0
  x <- pmax(x, eps)
  x <- x / sum(x)

  return(x)
}

softmax_to_simplex <- function(par_vec) {

  x <- c(par_vec, 0)
  x <- x - max(x)

  p <- exp(x)
  p <- p / sum(p)

  names(p) <- pi_names

  return(p)
}

simplex_to_softmax_start <- function(pi_vec, eps = 1e-12) {

  pi_vec <- project_simplex(pi_vec, eps = eps)

  log(pi_vec[1:7] / pi_vec[8])
}

calc_metrics <- function(est, true) {

  data.frame(
    MAE = mean(abs(est - true), na.rm = TRUE),
    RMSE = sqrt(mean((est - true)^2, na.rm = TRUE))
  )
}

make_candidate_row_no_true <- function(
  method,
  pi_vec,
  reference_pi,
  lambda_single = NA_real_,
  lambda_shared = NA_real_,
  lambda_zero = NA_real_,
  gamma = NA_real_,
  objective_value = NA_real_,
  convergence = NA_integer_,
  function_count = NA_real_,
  final_selected = FALSE,
  note = NA_character_
) {

  if (is.null(note) || length(note) == 0) {
    note <- NA_character_
  }

  if (is.null(objective_value) || length(objective_value) == 0) {
    objective_value <- NA_real_
  }

  if (is.null(convergence) || length(convergence) == 0) {
    convergence <- NA_integer_
  }

  if (is.null(function_count) || length(function_count) == 0) {
    function_count <- NA_real_
  }

  pi_vec <- project_simplex(pi_vec)
  names(pi_vec) <- pi_names

  reference_pi <- project_simplex(reference_pi)
  names(reference_pi) <- pi_names

  abs_change <- abs(pi_vec - reference_pi)
  rel_change <- (pi_vec - reference_pi) / pmax(abs(reference_pi), 1e-12)

  out <- data.frame(
    method = as.character(method)[1],
    lambda_single = as.numeric(lambda_single)[1],
    lambda_shared = as.numeric(lambda_shared)[1],
    lambda_zero = as.numeric(lambda_zero)[1],
    gamma = as.numeric(gamma)[1],
    objective_value = as.numeric(objective_value)[1],
    convergence = as.integer(convergence)[1],
    function_count = as.numeric(function_count)[1],
    total_abs_change_from_raw = sum(abs_change, na.rm = TRUE),
    max_abs_change_from_raw = max(abs_change, na.rm = TRUE),
    final_selected = final_selected,
    note = as.character(note)[1],
    stringsAsFactors = FALSE
  )

  for (nm in pi_names) {
    out[[nm]] <- as.numeric(pi_vec[nm])
  }

  for (nm in pi_names) {
    out[[paste0(nm, "_abs_change_from_raw")]] <- as.numeric(abs_change[nm])
  }

  for (nm in pi_names) {
    out[[paste0(nm, "_rel_change_from_raw")]] <- as.numeric(rel_change[nm])
  }

  out
}

make_candidate_row_comparison <- function(
  method,
  pi_vec,
  true_pi,
  lambda_single = NA_real_,
  lambda_shared = NA_real_,
  lambda_zero = NA_real_,
  gamma = NA_real_,
  objective_value = NA_real_,
  convergence = NA_integer_,
  function_count = NA_real_,
  final_selected = FALSE,
  note = NA_character_
) {

  if (is.null(note) || length(note) == 0) {
    note <- NA_character_
  }

  if (is.null(objective_value) || length(objective_value) == 0) {
    objective_value <- NA_real_
  }

  if (is.null(convergence) || length(convergence) == 0) {
    convergence <- NA_integer_
  }

  if (is.null(function_count) || length(function_count) == 0) {
    function_count <- NA_real_
  }

  pi_vec <- project_simplex(pi_vec)
  names(pi_vec) <- pi_names

  metric <- calc_metrics(pi_vec, true_pi)

  out <- data.frame(
    method = as.character(method)[1],
    lambda_single = as.numeric(lambda_single)[1],
    lambda_shared = as.numeric(lambda_shared)[1],
    lambda_zero = as.numeric(lambda_zero)[1],
    gamma = as.numeric(gamma)[1],
    objective_value = as.numeric(objective_value)[1],
    convergence = as.integer(convergence)[1],
    function_count = as.numeric(function_count)[1],
    final_selected = final_selected,
    MAE = as.numeric(metric$MAE)[1],
    RMSE = as.numeric(metric$RMSE)[1],
    note = as.character(note)[1],
    stringsAsFactors = FALSE
  )

  for (nm in pi_names) {
    out[[nm]] <- as.numeric(pi_vec[nm])
  }

  out
}

candidate_row_to_pi <- function(candidate_row) {

  pi_vec <- as.numeric(candidate_row[1, pi_names])
  names(pi_vec) <- pi_names
  project_simplex(pi_vec)
}

#============================================================
# LIKELIHOOD FUNCTIONS
#============================================================

neg_loglik_pi <- function(par_vec, soft_counts, eps = 1e-12) {

  pi_vec <- softmax_to_simplex(par_vec)

  pi_vec <- pmax(pi_vec, eps)
  pi_vec <- pi_vec / sum(pi_vec)

  value <- -sum(soft_counts * log(pi_vec), na.rm = TRUE)

  if (!is.finite(value)) {
    value <- .Machine$double.xmax
  }

  return(value)
}

neg_loglik_pi_directional <- function(
  par_vec,
  soft_counts,
  lambda_single,
  lambda_shared,
  eps = 1e-12
) {

  pi_vec <- softmax_to_simplex(par_vec)

  pi_vec <- pmax(pi_vec, eps)
  pi_vec <- pi_vec / sum(pi_vec)

  names(pi_vec) <- pi_names

  neg_ll <- -sum(soft_counts * log(pi_vec), na.rm = TRUE) /
    sum(soft_counts)

  penalty_single <-
    pi_vec["pi_1"] +
    pi_vec["pi_2"] +
    pi_vec["pi_3"]

  penalty_small_shared <- -log(pi_vec["pi_123"] + eps)

  value <-
    neg_ll +
    lambda_single * penalty_single +
    lambda_shared * penalty_small_shared

  if (!is.finite(value)) {
    value <- .Machine$double.xmax
  }

  return(value)
}

#============================================================
# OPTIMIZATION FUNCTIONS
#============================================================

optimize_pi_nelder_mead <- function(
  soft_counts,
  start_pi,
  maxit = 5000
) {

  start_pi <- project_simplex(start_pi)
  names(start_pi) <- pi_names

  par0 <- simplex_to_softmax_start(start_pi)

  opt <- optim(
    par = par0,
    fn = neg_loglik_pi,
    soft_counts = soft_counts,
    method = "Nelder-Mead",
    control = list(
      maxit = maxit,
      reltol = 1e-12
    )
  )

  pi_opt <- softmax_to_simplex(opt$par)
  pi_opt <- project_simplex(pi_opt)

  names(pi_opt) <- pi_names

  list(
    pi = pi_opt,
    value = opt$value,
    convergence = opt$convergence,
    counts = opt$counts,
    message = opt$message
  )
}

optimize_pi_nelder_mead_directional <- function(
  soft_counts,
  start_pi,
  lambda_single,
  lambda_shared,
  maxit = 5000
) {

  start_pi <- project_simplex(start_pi)
  names(start_pi) <- pi_names

  par0 <- simplex_to_softmax_start(start_pi)

  opt <- optim(
    par = par0,
    fn = neg_loglik_pi_directional,
    soft_counts = soft_counts,
    lambda_single = lambda_single,
    lambda_shared = lambda_shared,
    method = "Nelder-Mead",
    control = list(
      maxit = maxit,
      reltol = 1e-10
    )
  )

  pi_opt <- softmax_to_simplex(opt$par)
  pi_opt <- project_simplex(pi_opt)

  names(pi_opt) <- pi_names

  list(
    pi = pi_opt,
    value = opt$value,
    convergence = opt$convergence,
    counts = opt$counts,
    message = opt$message,
    lambda_single = lambda_single,
    lambda_shared = lambda_shared
  )
}

#============================================================
# DENSITY FILE FUNCTIONS
#============================================================

get_rho_index <- function(selected_rho, rho_vec) {

  idx <- match(
    round(selected_rho, 10),
    round(rho_vec, 10)
  )

  if (is.na(idx)) {
    stop("selected_rho not found in rho_vec: ", selected_rho)
  }

  idx
}

build_density_file <- function(chr, sim, tauuse, density_dir) {

  file.path(
    density_dir,
    paste0(
      penalty,
      "_chr", chr,
      "_3pop_",
      "_warmStart", warmStart,
      "_sim", sim,
      "_Zscale", Zscale,
      "_singleStart", singleStart,
      "_tauuse", tauuse,
      "_DensityU.RData"
    )
  )
}

density_row_to_soft_counts <- function(density_row) {

  if (length(density_row) %% 8 != 0) {
    stop("Length of density_row is not divisible by 8")
  }

  P <- length(density_row) / 8

  u_mat <- matrix(
    density_row,
    nrow = P,
    byrow = TRUE
  )

  row_sums <- rowSums(u_mat)
  row_sums[!is.finite(row_sums) | row_sums <= 0] <- 1

  post_prob <- u_mat / row_sums
  post_prob[!is.finite(post_prob)] <- 0

  colnames(post_prob) <- pi_names

  list(
    soft_counts = colSums(post_prob),
    P = P
  )
}

#============================================================
# CORRECTION FUNCTIONS
#============================================================

correct_pi123_v2 <- function(
  pi_vec,
  lambda_single,
  lambda_zero
) {

  pi_vec <- project_simplex(pi_vec)
  names(pi_vec) <- pi_names

  pi_new <- pi_vec

  single_names <- c("pi_1", "pi_2", "pi_3")

  single_mass <- sum(pi_vec[single_names], na.rm = TRUE)
  zero_mass <- pi_vec["pi_0"]

  move_from_single <- lambda_single * single_mass
  move_from_zero <- lambda_zero * zero_mass

  # Population-aware correction:
  # pi_1 = AFR, pi_2 = EAS, pi_3 = EUR
  # Shrink each ancestry-specific component proportionally.
  # This avoids artificially forcing the smallest ancestry component
  # especially EAS/pi_2 to zero.
  pi_new[single_names] <- pi_vec[single_names] * (1 - lambda_single)

  pi_new["pi_0"] <- pi_vec["pi_0"] - move_from_zero

  pi_new["pi_123"] <- pi_vec["pi_123"] +
    move_from_single +
    move_from_zero

  pi_new <- project_simplex(pi_new)
  names(pi_new) <- names(pi_vec)

  return(pi_new)
}


run_correction_grid_no_true <- function(
  pi_vec,
  correction_grid
) {

  grid_rows <- list()

  for (ii in seq_len(nrow(correction_grid))) {

    lambda_single <- correction_grid$lambda_single[ii]
    lambda_zero <- correction_grid$lambda_zero[ii]

    pi_corr <- correct_pi123_v2(
      pi_vec = pi_vec,
      lambda_single = lambda_single,
      lambda_zero = lambda_zero
    )

    distortion <- sum((pi_corr - pi_vec)^2, na.rm = TRUE)

    row <- make_candidate_row_no_true(
      method = "corrected_pi123_grid_diagnostic_no_true",
      pi_vec = pi_corr,
      reference_pi = pi_vec,
      lambda_single = lambda_single,
      lambda_zero = lambda_zero,
      objective_value = distortion,
      convergence = 0,
      function_count = NA_real_,
      final_selected = FALSE,
      note = "Diagnostic only; final lambda fixed from simulation calibration."
    )

    row$distortion_from_raw_pi_final <- distortion
    row$pi123_change_from_raw <- as.numeric(pi_corr["pi_123"] - pi_vec["pi_123"])

    grid_rows[[length(grid_rows) + 1]] <- row
  }

  do.call(rbind, grid_rows)
}

run_directional_lambda_grid_no_true <- function(
  soft_counts,
  start_pi,
  directional_lambda_grid
) {

  grid_rows <- list()

  for (ii in seq_len(nrow(directional_lambda_grid))) {

    lambda_single <- directional_lambda_grid$lambda_single[ii]
    lambda_shared <- directional_lambda_grid$lambda_shared[ii]

    opt <- optimize_pi_nelder_mead_directional(
      soft_counts = soft_counts,
      start_pi = start_pi,
      lambda_single = lambda_single,
      lambda_shared = lambda_shared
    )

    grid_rows[[length(grid_rows) + 1]] <- make_candidate_row_no_true(
      method = "directional_nelder_mead_grid_diagnostic_no_true",
      pi_vec = opt$pi,
      reference_pi = start_pi,
      lambda_single = lambda_single,
      lambda_shared = lambda_shared,
      objective_value = opt$value,
      convergence = opt$convergence,
      function_count = as.numeric(opt$counts["function"]),
      final_selected = FALSE,
      note = ifelse(
        is.null(opt$message) || length(opt$message) == 0,
        NA_character_,
        opt$message
      )
    )
  }

  do.call(rbind, grid_rows)
}

run_unsupervised_lambda_score_grid <- function(
  pi_vec,
  correction_grid,
  gamma_grid
) {

  grid_rows <- list()

  for (ii in seq_len(nrow(correction_grid))) {

    lambda_single <- correction_grid$lambda_single[ii]
    lambda_zero <- correction_grid$lambda_zero[ii]

    pi_corr <- correct_pi123_v2(
      pi_vec = pi_vec,
      lambda_single = lambda_single,
      lambda_zero = lambda_zero
    )

    distortion <- sum((pi_corr - pi_vec)^2, na.rm = TRUE)
    gain_123 <- pi_corr["pi_123"]

    for (gamma in gamma_grid) {

      score <- gain_123 - gamma * distortion

      grid_rows[[length(grid_rows) + 1]] <- data.frame(
        lambda_single = lambda_single,
        lambda_zero = lambda_zero,
        gamma = gamma,
        pi_123 = as.numeric(gain_123),
        distortion = distortion,
        score = as.numeric(score),
        note = "Diagnostic score only; not used for final selection.",
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, grid_rows)
}

#============================================================
# SIMULATION COMPARISON FUNCTIONS
#============================================================
# These functions use true_pi only for comparison. They do not choose the final estimate.

make_comparison_table <- function(
  pi_final,
  pi_optimized_nm,
  pi_final_optimized,
  pi_directional_calibrated,
  pi_es_alpha,
  true_pi
) {

  rows <- list()

  rows[[length(rows) + 1]] <- make_candidate_row_comparison(
    method = "pi_final_raw",
    pi_vec = pi_final,
    true_pi = true_pi,
    note = "Comparison only."
  )

  rows[[length(rows) + 1]] <- make_candidate_row_comparison(
    method = "basic_nelder_mead_check",
    pi_vec = pi_optimized_nm,
    true_pi = true_pi,
    note = "Comparison only."
  )

  rows[[length(rows) + 1]] <- make_candidate_row_comparison(
    method = "final_calibrated_correct_pi123_v2_no_true",
    pi_vec = pi_final_optimized,
    true_pi = true_pi,
    lambda_single = lambda_single_calibrated,
    lambda_zero = lambda_zero_calibrated,
    final_selected = TRUE,
    note = "Final estimate already computed without true_pi; true_pi used here only for comparison."
  )

  rows[[length(rows) + 1]] <- make_candidate_row_comparison(
    method = "directional_calibrated_comparison",
    pi_vec = pi_directional_calibrated,
    true_pi = true_pi,
    lambda_single = lambda_single_directional_calibrated,
    lambda_shared = lambda_shared_directional_calibrated,
    note = "Comparison only."
  )

  if (all(is.finite(pi_es_alpha))) {
    rows[[length(rows) + 1]] <- make_candidate_row_comparison(
      method = "pi_es_alpha",
      pi_vec = pi_es_alpha,
      true_pi = true_pi,
      note = "Comparison only."
    )
  }

  do.call(rbind, rows)
}

#============================================================
# READ SUMMARY FILES
#============================================================

inputs_list <- lapply(summary_file, process_one_summary)
inputs_list <- Filter(Negate(is.null), inputs_list)

if (length(inputs_list) == 0) {
  stop("All summary files were empty or invalid in: ", sc_dir)
}

summary_df <- do.call(
  rbind,
  lapply(inputs_list, as.data.frame)
)

rownames(summary_df) <- NULL

cat("Scenario:", sc_name, "\n")
cat("Scenario index:", scenarioIndex, "\n")
cat("Usable summary files:", nrow(summary_df), "\n")

print(summary_df)

#============================================================
# MAIN LOOP
#============================================================

results_all <- list()
metrics_all <- list()
candidate_summary_all <- list()
correction_grid_no_true_all <- list()
directional_grid_no_true_all <- list()
unsupervised_score_grid_all <- list()
comparison_summary_all <- list()

for (i in seq_len(nrow(summary_df))) {

  sim <- summary_df$sim[i]
  selected_roww <- summary_df$selected_roww[i]
  selected_rho <- summary_df$selected_rho[i]
  tauuse <- summary_df$tauuse[i]
  largest_logsum <- summary_df$Largest_Sum_Of_Logs[i]

  cat("\n============================================================\n")
  cat(
    "Processing sim =", sim,
    "| selected_row =", selected_roww,
    "| selected_rho =", selected_rho,
    "| tauuse =", tauuse,
    "\n"
  )
  cat("============================================================\n")

  #============================================================
  # SAFE ES-ALPHA LOADING
  #============================================================

  es_alpha_file <- file.path(
    sc_dir,
    paste0(
      "es_alpha_best_sim", sim,
      "_scenario", scenarioIndex,
      "_row", selected_roww,
      "_tau", tauuse,
      ".RData"
    )
  )

  if (file.exists(es_alpha_file)) {

    e_es <- new.env(parent = emptyenv())

    ok_es <- tryCatch(
      {
        load(es_alpha_file, envir = e_es)
        TRUE
      },
      error = function(err) {
        warning("Could not load es_alpha file: ", es_alpha_file)
        FALSE
      }
    )

    if (ok_es && exists("selected_row", envir = e_es, inherits = FALSE)) {

      selected_row_es <- get("selected_row", envir = e_es)

      if (length(selected_row_es) == 8) {
        pi_es_alpha <- project_simplex(selected_row_es)
      } else {
        warning("selected_row missing or invalid inside: ", es_alpha_file)
        pi_es_alpha <- rep(NA_real_, 8)
      }

    } else {
      warning("selected_row missing inside: ", es_alpha_file)
      pi_es_alpha <- rep(NA_real_, 8)
    }

    rm(e_es)

  } else {

    warning("Missing es_alpha for sim ", sim)
    pi_es_alpha <- rep(NA_real_, 8)
  }

  names(pi_es_alpha) <- pi_names

  #============================================================
  # POSTERIOR COUNTS
  #============================================================

  rho_idx <- get_rho_index(selected_rho, rho_vec)

  pooled_counts <- rep(0, 8)
  names(pooled_counts) <- pi_names

  total_snps_used <- 0
  valid_chr <- 0

  for (chr in 1:22) {

    print(chr)

    f <- build_density_file(
      chr = chr,
      sim = sim,
      tauuse = tauuse,
      density_dir = density_dir
    )

    if (!file.exists(f)) next
    if (file.info(f)$size == 0) next

    e <- new.env(parent = emptyenv())

    ok <- tryCatch(
      {
        load(f, envir = e)
        TRUE
      },
      error = function(err) {
        message("Skipping corrupted file: ", f)
        FALSE
      }
    )

    if (!ok) {
      rm(e)
      next
    }

    if (!exists("Allrho_DensityU_list", envir = e, inherits = FALSE)) {
      rm(e)
      next
    }

    dens_list <- get("Allrho_DensityU_list", envir = e)

    if (length(dens_list) < rho_idx) {
      rm(e)
      next
    }

    dens_mat <- dens_list[[rho_idx]]

    if (is.null(dens_mat)) {
      rm(e)
      next
    }

    if (selected_roww < 1 || selected_roww > nrow(dens_mat)) {
      stop(
        "selected_row out of bounds for chr ",
        chr,
        " in sim ",
        sim
      )
    }

    one <- density_row_to_soft_counts(
      dens_mat[selected_roww, ]
    )

    pooled_counts <- pooled_counts + one$soft_counts
    total_snps_used <- total_snps_used + one$P
    valid_chr <- valid_chr + 1

    rm(e)
    gc()
  }

  if (sum(pooled_counts) <= 0) {
    warning("No valid posterior counts for sim ", sim)
    next
  }

  pi_final <- project_simplex(pooled_counts)
  names(pi_final) <- pi_names

  #============================================================
  # BASIC NELDER-MEAD OPTIMIZATION FOR PI
  #============================================================

  nm_from_final <- optimize_pi_nelder_mead(
    soft_counts = pooled_counts,
    start_pi = pi_final
  )

  if (all(is.finite(pi_es_alpha))) {
    nm_from_es <- optimize_pi_nelder_mead(
      soft_counts = pooled_counts,
      start_pi = pi_es_alpha
    )
  } else {
    nm_from_es <- NULL
  }

  if (!is.null(nm_from_es) && nm_from_es$value < nm_from_final$value) {

    pi_optimized_nm <- nm_from_es$pi
    nm_start_used <- "pi_es_alpha"
    nm_value <- nm_from_es$value
    nm_convergence <- nm_from_es$convergence
    nm_function_count <- nm_from_es$counts["function"]

  } else {

    pi_optimized_nm <- nm_from_final$pi
    nm_start_used <- "pi_final"
    nm_value <- nm_from_final$value
    nm_convergence <- nm_from_final$convergence
    nm_function_count <- nm_from_final$counts["function"]
  }

  names(pi_optimized_nm) <- pi_names

  #============================================================
  # DIAGNOSTIC GRIDS WITHOUT TRUE VALUES
  #============================================================
  # These grids show sensitivity, but they do not select the final estimate.

  correction_grid_df <- run_correction_grid_no_true(
    pi_vec = pi_final,
    correction_grid = correction_grid
  )

  directional_grid_df <- run_directional_lambda_grid_no_true(
    soft_counts = pooled_counts,
    start_pi = pi_final,
    directional_lambda_grid = directional_lambda_grid
  )

  unsupervised_score_grid_df <- run_unsupervised_lambda_score_grid(
    pi_vec = pi_final,
    correction_grid = correction_grid,
    gamma_grid = gamma_grid
  )

  correction_grid_df$sim <- sim
  correction_grid_df$scenarioIndex <- scenarioIndex
  correction_grid_df$scenarioName <- sc_name

  directional_grid_df$sim <- sim
  directional_grid_df$scenarioIndex <- scenarioIndex
  directional_grid_df$scenarioName <- sc_name

  unsupervised_score_grid_df$sim <- sim
  unsupervised_score_grid_df$scenarioIndex <- scenarioIndex
  unsupervised_score_grid_df$scenarioName <- sc_name

  #============================================================
  # FINAL ESTIMATE WITHOUT USING TRUE VALUES
  #============================================================
  # This is the final practical estimator for real data.

  pi_final_optimized <- correct_pi123_v2(
    pi_vec = pi_final,
    lambda_single = lambda_single_calibrated,
    lambda_zero = lambda_zero_calibrated
  )

  pi_corrected_123 <- pi_final_optimized

  # Optional directional comparison using fixed calibrated values.
  # Not final unless separately validated.

  directional_calibrated_out <- optimize_pi_nelder_mead_directional(
    soft_counts = pooled_counts,
    start_pi = pi_final,
    lambda_single = lambda_single_directional_calibrated,
    lambda_shared = lambda_shared_directional_calibrated
  )

  pi_optimized_directional_nm <- directional_calibrated_out$pi
  names(pi_optimized_directional_nm) <- pi_names

  #============================================================
  # CANDIDATE SUMMARY WITHOUT TRUE VALUES
  #============================================================

  candidate_grid_df <- rbind(
    make_candidate_row_no_true(
      method = "pi_final_raw",
      pi_vec = pi_final,
      reference_pi = pi_final,
      objective_value = NA_real_,
      convergence = 0,
      final_selected = FALSE,
      note = "Raw posterior soft-count estimate."
    ),
    make_candidate_row_no_true(
      method = "basic_nelder_mead_check",
      pi_vec = pi_optimized_nm,
      reference_pi = pi_final,
      objective_value = nm_value,
      convergence = nm_convergence,
      function_count = as.numeric(nm_function_count),
      final_selected = FALSE,
      note = nm_start_used
    ),
    make_candidate_row_no_true(
      method = "final_calibrated_correct_pi123_v2_no_true",
      pi_vec = pi_final_optimized,
      reference_pi = pi_final,
      lambda_single = lambda_single_calibrated,
      lambda_zero = lambda_zero_calibrated,
      objective_value = sum((pi_final_optimized - pi_final)^2, na.rm = TRUE),
      convergence = 0,
      final_selected = TRUE,
      note = "Final estimate; lambdas fixed from simulation calibration, not selected by true_pi."
    ),
    make_candidate_row_no_true(
      method = "directional_calibrated_comparison",
      pi_vec = pi_optimized_directional_nm,
      reference_pi = pi_final,
      lambda_single = lambda_single_directional_calibrated,
      lambda_shared = lambda_shared_directional_calibrated,
      objective_value = directional_calibrated_out$value,
      convergence = directional_calibrated_out$convergence,
      function_count = as.numeric(directional_calibrated_out$counts["function"]),
      final_selected = FALSE,
      note = "Comparison only unless validated by simulation calibration."
    )
  )

  if (all(is.finite(pi_es_alpha))) {
    candidate_grid_df <- rbind(
      candidate_grid_df,
      make_candidate_row_no_true(
        method = "pi_es_alpha",
        pi_vec = pi_es_alpha,
        reference_pi = pi_final,
        objective_value = NA_real_,
        convergence = 0,
        final_selected = FALSE,
        note = "External ES-alpha comparison."
      )
    )
  }

  candidate_grid_df$sim <- sim
  candidate_grid_df$scenarioIndex <- scenarioIndex
  candidate_grid_df$scenarioName <- sc_name

  #============================================================
  # SIMULATION COMPARISON ONLY
  #============================================================
  # true_pi is used only below this line.
  # It does not change pi_final_optimized.

  if (use_true_pi_for_comparison) {

    metric_final <- calc_metrics(pi_final, true_pi)
    metric_es <- calc_metrics(pi_es_alpha, true_pi)
    metric_corrected_123 <- calc_metrics(pi_corrected_123, true_pi)
    metric_optimized_nm <- calc_metrics(pi_optimized_nm, true_pi)
    metric_optimized_directional_nm <- calc_metrics(
      pi_optimized_directional_nm,
      true_pi
    )
    metric_final_optimized <- calc_metrics(
      pi_final_optimized,
      true_pi
    )

    comparison_summary_df <- make_comparison_table(
      pi_final = pi_final,
      pi_optimized_nm = pi_optimized_nm,
      pi_final_optimized = pi_final_optimized,
      pi_directional_calibrated = pi_optimized_directional_nm,
      pi_es_alpha = pi_es_alpha,
      true_pi = true_pi
    )

    comparison_summary_df$sim <- sim
    comparison_summary_df$scenarioIndex <- scenarioIndex
    comparison_summary_df$scenarioName <- sc_name

    true_out <- as.numeric(true_pi)
    bias_final <- as.numeric(pi_final - true_pi)
    bias_nm <- as.numeric(pi_optimized_nm - true_pi)
    bias_directional <- as.numeric(pi_optimized_directional_nm - true_pi)
    bias_corrected <- as.numeric(pi_corrected_123 - true_pi)
    bias_final_optimized <- as.numeric(pi_final_optimized - true_pi)
    bias_es <- as.numeric(pi_es_alpha - true_pi)

  } else {

    metric_final <- data.frame(MAE = NA_real_, RMSE = NA_real_)
    metric_es <- data.frame(MAE = NA_real_, RMSE = NA_real_)
    metric_corrected_123 <- data.frame(MAE = NA_real_, RMSE = NA_real_)
    metric_optimized_nm <- data.frame(MAE = NA_real_, RMSE = NA_real_)
    metric_optimized_directional_nm <- data.frame(MAE = NA_real_, RMSE = NA_real_)
    metric_final_optimized <- data.frame(MAE = NA_real_, RMSE = NA_real_)

    comparison_summary_df <- data.frame()

    true_out <- rep(NA_real_, length(pi_names))
    bias_final <- rep(NA_real_, length(pi_names))
    bias_nm <- rep(NA_real_, length(pi_names))
    bias_directional <- rep(NA_real_, length(pi_names))
    bias_corrected <- rep(NA_real_, length(pi_names))
    bias_final_optimized <- rep(NA_real_, length(pi_names))
    bias_es <- rep(NA_real_, length(pi_names))
  }

  #============================================================
  # RESULT TABLE
  #============================================================

  sim_result <- data.frame(
    sim = sim,
    scenarioIndex = scenarioIndex,
    scenarioName = sc_name,
    selected_row = selected_roww,
    selected_rho = selected_rho,
    tauuse = tauuse,
    Largest_Sum_Of_Logs = largest_logsum,
    Parameter = pi_names,
    True_Comparison_Only = true_out,
    Pi_Final = as.numeric(pi_final),
    Pi_Optimized_NelderMead = as.numeric(pi_optimized_nm),
    Pi_Directional_Calibrated_Comparison =
      as.numeric(pi_optimized_directional_nm),
    Pi_Best_Corrected_123_NoTrue = as.numeric(pi_corrected_123),
    Pi_Final_Optimized_NoTrue = as.numeric(pi_final_optimized),
    Pi_es_alpha = as.numeric(pi_es_alpha),
    Bias_Final_Comparison_Only = bias_final,
    Bias_Optimized_NelderMead_Comparison_Only = bias_nm,
    Bias_Directional_Calibrated_Comparison_Only = bias_directional,
    Bias_Best_Corrected_123_NoTrue_Comparison_Only = bias_corrected,
    Bias_Final_Optimized_NoTrue_Comparison_Only = bias_final_optimized,
    Bias_es_alpha_Comparison_Only = bias_es,
    AbsBias_Final_Comparison_Only = abs(bias_final),
    AbsBias_Optimized_NelderMead_Comparison_Only = abs(bias_nm),
    AbsBias_Directional_Calibrated_Comparison_Only = abs(bias_directional),
    AbsBias_Best_Corrected_123_NoTrue_Comparison_Only = abs(bias_corrected),
    AbsBias_Final_Optimized_NoTrue_Comparison_Only = abs(bias_final_optimized),
    AbsBias_es_alpha_Comparison_Only = abs(bias_es),
    Final_Estimator_Note = "Pi_Final_Optimized_NoTrue is computed using fixed calibrated lambdas, not true_pi.",
    stringsAsFactors = FALSE
  )

  sim_metrics <- data.frame(
    sim = sim,
    scenarioIndex = scenarioIndex,
    scenarioName = sc_name,
    selected_row = selected_roww,
    selected_rho = selected_rho,
    tauuse = tauuse,
    Largest_Sum_Of_Logs = largest_logsum,
    total_snps_used = total_snps_used,
    valid_chr = valid_chr,
    NelderMead_Start = nm_start_used,
    NelderMead_NegLogLik = nm_value,
    NelderMead_Convergence = nm_convergence,
    NelderMead_Function_Count = as.numeric(nm_function_count),
    Final_Optimized_Method = "final_calibrated_correct_pi123_v2_no_true",
    Final_Optimized_Lambda_Single = lambda_single_calibrated,
    Final_Optimized_Lambda_Shared = NA_real_,
    Final_Optimized_Lambda_Zero = lambda_zero_calibrated,
    Final_Optimized_Gamma = NA_real_,
    Final_Optimized_Objective =
      sum((pi_final_optimized - pi_final)^2, na.rm = TRUE),
    Final_Optimized_Convergence = 0,
    Final_Optimized_Function_Count = NA_real_,
    Final_Parameter_Choice_Note =
      "Final lambdas fixed from aggregate simulation calibration across 5 scenarios; true_pi not used for selection.",
    Directional_Calibrated_Lambda_Single = lambda_single_directional_calibrated,
    Directional_Calibrated_Lambda_Shared = lambda_shared_directional_calibrated,
    Directional_Calibrated_Objective = directional_calibrated_out$value,
    Directional_Calibrated_Convergence = directional_calibrated_out$convergence,
    Directional_Calibrated_Function_Count =
      as.numeric(directional_calibrated_out$counts["function"]),
    Total_Abs_Change_Optimized_vs_Raw =
      sum(abs(pi_final_optimized - pi_final), na.rm = TRUE),
    Max_Abs_Change_Optimized_vs_Raw =
      max(abs(pi_final_optimized - pi_final), na.rm = TRUE),
    Pi123_Raw = as.numeric(pi_final["pi_123"]),
    Pi123_Optimized_NoTrue = as.numeric(pi_final_optimized["pi_123"]),
    MAE_Final_Comparison_Only = metric_final$MAE,
    RMSE_Final_Comparison_Only = metric_final$RMSE,
    MAE_Optimized_NelderMead_Comparison_Only = metric_optimized_nm$MAE,
    RMSE_Optimized_NelderMead_Comparison_Only = metric_optimized_nm$RMSE,
    MAE_Directional_Calibrated_Comparison_Only =
      metric_optimized_directional_nm$MAE,
    RMSE_Directional_Calibrated_Comparison_Only =
      metric_optimized_directional_nm$RMSE,
    MAE_Best_Corrected_123_NoTrue_Comparison_Only =
      metric_corrected_123$MAE,
    RMSE_Best_Corrected_123_NoTrue_Comparison_Only =
      metric_corrected_123$RMSE,
    MAE_Final_Optimized_NoTrue_Comparison_Only =
      metric_final_optimized$MAE,
    RMSE_Final_Optimized_NoTrue_Comparison_Only =
      metric_final_optimized$RMSE,
    MAE_es_alpha_Comparison_Only = metric_es$MAE,
    RMSE_es_alpha_Comparison_Only = metric_es$RMSE,
    stringsAsFactors = FALSE
  )

  #============================================================
  # PRINT RESULTS
  #============================================================

  cat("\n==== BASIC NELDER-MEAD CHECK FOR SIM", sim, "====\n")
  cat("Start used:", nm_start_used, "\n")
  cat("Negative log-likelihood:", nm_value, "\n")
  cat("Convergence code:", nm_convergence, "\n")

  cat("\n==== FINAL ESTIMATE WITHOUT TRUE_PI FOR SIM", sim, "====\n")
  print(
    candidate_grid_df[
      candidate_grid_df$final_selected == TRUE,
      c(
        "method",
        "lambda_single",
        "lambda_zero",
        "objective_value",
        "total_abs_change_from_raw",
        "max_abs_change_from_raw",
        "final_selected",
        pi_names
      ),
      drop = FALSE
    ]
  )

  if (use_true_pi_for_comparison) {
    cat("\n==== SIMULATION COMPARISON ONLY FOR SIM", sim, "====\n")
    print(
      comparison_summary_df[
        ,
        c(
          "method",
          "lambda_single",
          "lambda_shared",
          "lambda_zero",
          "final_selected",
          "MAE",
          "RMSE",
          "pi_123"
        ),
        drop = FALSE
      ]
    )
  }

  cat("\n==== FINAL PI COMPARISON TABLE FOR SIM", sim, "====\n")
  print(
    sim_result[
      ,
      c(
        "Parameter",
        "True_Comparison_Only",
        "Pi_Final",
        "Pi_Optimized_NelderMead",
        "Pi_Directional_Calibrated_Comparison",
        "Pi_Best_Corrected_123_NoTrue",
        "Pi_Final_Optimized_NoTrue",
        "Pi_es_alpha",
        "Bias_Final_Optimized_NoTrue_Comparison_Only",
        "Bias_es_alpha_Comparison_Only"
      )
    ]
  )

  cat("\n==== ERROR SUMMARY / COMPARISON ONLY FOR SIM", sim, "====\n")
  print(sim_metrics)

  results_all[[length(results_all) + 1]] <- sim_result
  metrics_all[[length(metrics_all) + 1]] <- sim_metrics
  candidate_summary_all[[length(candidate_summary_all) + 1]] <- candidate_grid_df
  correction_grid_no_true_all[[length(correction_grid_no_true_all) + 1]] <-
    correction_grid_df
  directional_grid_no_true_all[[length(directional_grid_no_true_all) + 1]] <-
    directional_grid_df
  unsupervised_score_grid_all[[length(unsupervised_score_grid_all) + 1]] <-
    unsupervised_score_grid_df

  if (use_true_pi_for_comparison) {
    comparison_summary_all[[length(comparison_summary_all) + 1]] <-
      comparison_summary_df
  }
}

#============================================================
# SAVE FINAL OUTPUTS
#============================================================

if (length(results_all) == 0) {
  stop("No simulation produced valid results.")
}

final_results <- do.call(rbind, results_all)
final_metrics <- do.call(rbind, metrics_all)
final_candidate_summary <- do.call(rbind, candidate_summary_all)
final_correction_grid_no_true <- do.call(rbind, correction_grid_no_true_all)
final_directional_grid_no_true <- do.call(rbind, directional_grid_no_true_all)
final_unsupervised_score_grid <- do.call(rbind, unsupervised_score_grid_all)

if (length(comparison_summary_all) > 0) {
  final_comparison_summary <- do.call(rbind, comparison_summary_all)
} else {
  final_comparison_summary <- data.frame()
}

out_results <- paste0(
  sc_dir,
  "/Final_pi_notrue_results.csv"
)

out_metrics <- paste0(
  sc_dir,
  "/Final_pi_notrue_metrics_with_comparison.csv"
)

out_candidate_summary <- paste0(
  sc_dir,
  "/Final_pi_candidate_summary_NO_TRUE_USED_FOR_SELECTION.csv"
)

out_correction_grid <- paste0(
  sc_dir,
  "/Final_pi_correction_grid_diagnostic_NO_TRUE.csv"
)

out_directional_grid <- paste0(
  sc_dir,
  "/Final_pi_directional_grid_diagnostic_NO_TRUE.csv"
)

out_unsupervised_score_grid <- paste0(
  sc_dir,
  "/Final_pi_unsupervised_lambda_score_grid_NO_TRUE.csv"
)

out_comparison_summary <- paste0(
  sc_dir,
  "/Final_pi_TRUE_COMPARISON_ONLY_summary.csv"
)

write.csv(final_results, out_results, row.names = FALSE)
write.csv(final_metrics, out_metrics, row.names = FALSE)
write.csv(final_candidate_summary, out_candidate_summary, row.names = FALSE)
write.csv(final_correction_grid_no_true, out_correction_grid, row.names = FALSE)
write.csv(final_directional_grid_no_true, out_directional_grid, row.names = FALSE)
write.csv(final_unsupervised_score_grid, out_unsupervised_score_grid, row.names = FALSE)

if (nrow(final_comparison_summary) > 0) {
  write.csv(final_comparison_summary, out_comparison_summary, row.names = FALSE)
}

cat("\nSaved files:\n")
cat(out_results, "\n")
cat(out_metrics, "\n")
cat(out_candidate_summary, "\n")
cat(out_correction_grid, "\n")
cat(out_directional_grid, "\n")
cat(out_unsupervised_score_grid, "\n")

if (nrow(final_comparison_summary) > 0) {
  cat(out_comparison_summary, "\n")
}

cat("\nIMPORTANT:\n")
cat("Pi_Final_Optimized_NoTrue is the final estimator.\n")
cat("true_pi is used only for comparison columns ending with _Comparison_Only.\n")
cat("No RMSE/MAE/Bias value is used to choose the final estimate.\n")
