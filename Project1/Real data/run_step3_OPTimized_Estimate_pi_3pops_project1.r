rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

#============================================================
# REAL-DATA-SAFE PI ESTIMATION SCRIPT
#============================================================
# Goal:
#   1. Estimate pi from posterior soft counts using real DensityU files.
#   2. Compute ES-alpha from real beta matrices across all chromosomes.
#   3. Compute a simple simulation-calibrated correction as sensitivity.
#   4. Compute the directional simulation-calibrated estimate as FINAL.
#   5. NEVER use true_pi for real data.
#
# Important:
#   - Pi_Final_Optimized_NoTrue is the practical FINAL estimate.
#   - In this corrected script, Pi_Final_Optimized_NoTrue is the
#     directional calibrated estimate.
#   - Pi_Best_Corrected_123_NoTrue is kept as a secondary sensitivity estimate.
#   - Calibrated lambdas are fixed before applying this script to real data.
#   - This script is for real data, so use_true_pi_for_comparison <- FALSE.
#============================================================

#============================================================
# INPUTS: REAL DATA
#============================================================

penalty <- "RealmixLOG"
savename <- "new_densityU_3pops_project1"

warmStart <- 1
Zscale <- 1
singleStart <- 1
K <- 3

popvec <- c("AFR", "EAS", "EUR")

dirOutput <- "/home/yatah3/projects/def-thchlava/yatah3/real_data_backup/real-data/project1/"

output_sub_folder <- "/lustre06/project/6005709/yatah3/real_data_backup/real-data/project1/Trans2-new-3-out/"

output_sub_folder_uDensity <- paste0(
  "/lustre06/project/6005709/yatah3/real_data_backup/real-data/project1/",
  savename,
  "/"
)

functionsfolder <- "/lustre06/project/6005709/yatah3/real_data_backup/real-data/Rfunction/"

dir.create(dirOutput, recursive = TRUE, showWarnings = FALSE)

#============================================================
# LOAD REQUIRED FUNCTIONS
#============================================================

source(paste0(functionsfolder, "AllPRS_Rfunctions.r"))
source(paste0(functionsfolder, "PRS_utility.r"))
source(paste0(functionsfolder, "Iterative_Rfunctions.r"))
source(paste0(functionsfolder, "PlinkLD_transform.R"))

#============================================================
# LOAD TAU FILE AND SELECT REAL-DATA TARGET TAU/RHO/ROW
#============================================================

Taufile <- paste0(
  dirOutput,
  popvec[1], "_", popvec[2], "-", popvec[3],
  "Tau_info_v2.RData"
)

if (!file.exists(Taufile)) {
  stop("Tau file not found: ", Taufile)
}

load(Taufile)

if (!exists("AbsTauvec")) {
  stop("AbsTauvec not found in Tau file: ", Taufile)
}

tau_target <- 7.5e-05
tauuse <- AbsTauvec[which.min(abs(AbsTauvec - tau_target))]

rho_vec <- c(seq(0, 0.9, 0.1), 0.95)
rho_target <- 0.95
selected_rho <- rho_target

selected_roww <- 68

cat("============================================================\n")
cat("REAL DATA SETTINGS\n")
cat("Population order:", paste(popvec, collapse = ", "), "\n")
cat("Selected row:", selected_roww, "\n")
cat("Selected rho:", selected_rho, "\n")
cat("Selected tau:", tauuse, "\n")
cat("Density directory:", output_sub_folder_uDensity, "\n")
cat("Beta directory:", output_sub_folder, "\n")
cat("Output directory:", dirOutput, "\n")
cat("============================================================\n")

#============================================================
# SETTINGS
#============================================================

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
# TRUE PI: REAL DATA
#============================================================
# For real data, true pi is unknown.
# It is never used in this script.

use_true_pi_for_comparison <- FALSE

true_pi <- rep(NA_real_, length(pi_names))
names(true_pi) <- pi_names

#============================================================
# CALIBRATED PARAMETERS FOR REAL-DATA-SAFE ESTIMATES
#============================================================
# These values are fixed from real-data-like simulation calibration.
# true_pi is not used anywhere in this real-data script.

# Secondary sensitivity estimate:
# Simple corrected pi123 rule.
lambda_single_calibrated <- 0.50
lambda_zero_calibrated   <- 0.00

# Primary final estimate:
# Directional simulation-calibrated method.
lambda_single_directional_calibrated <- 5
lambda_shared_directional_calibrated <- 0.025

#============================================================
# DIAGNOSTIC PARAMETER GRIDS
#============================================================
# Saved for sensitivity analysis only.
# Not used to select the final estimate in real data.

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

project_simplex <- function(x, eps = 1e-12) {

  x <- as.numeric(x)
  names(x) <- pi_names[seq_along(x)]

  x[!is.finite(x)] <- 0
  x <- pmax(x, eps)
  x <- x / sum(x)

  names(x) <- pi_names

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

  rownames(out) <- NULL

  return(out)
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

  return(idx)
}

find_density_file <- function(chr, tauuse, density_dir) {

  candidate_files <- c(
    file.path(
      density_dir,
      paste0(
        penalty,
        "_chr", chr,
        "_3pop_",
        "_warmStart", warmStart,
        "_Zscale", Zscale,
        "_singleStart", singleStart,
        "_tauuse", tauuse,
        "_DensityU.RData"
      )
    ),
    file.path(
      density_dir,
      paste0(
        penalty,
        "_chr", chr,
        "_3pops_",
        "_warmStart", warmStart,
        "_Zscale", Zscale,
        "_singleStart", singleStart,
        "_tauuse", tauuse,
        "_DensityU.RData"
      )
    ),
    file.path(
      density_dir,
      paste0(
        penalty,
        "_chr", chr,
        "_3traits",
        "warmStart", warmStart,
        "Zscale", Zscale,
        "tauuse", tauuse,
        "singleStart", singleStart,
        "_DensityU.RData"
      )
    )
  )

  existing <- candidate_files[file.exists(candidate_files)]

  if (length(existing) > 0) {
    return(existing[1])
  }

  pattern_chr <- paste0("chr", chr)

  all_files <- list.files(
    density_dir,
    pattern = "\\.RData$",
    full.names = TRUE
  )

  if (length(all_files) == 0) {
    return(NA_character_)
  }

  matched <- all_files[
    grepl(pattern_chr, basename(all_files)) &
      grepl("DensityU", basename(all_files)) &
      grepl(as.character(tauuse), basename(all_files), fixed = TRUE)
  ]

  if (length(matched) == 0) {
    matched <- all_files[
      grepl(pattern_chr, basename(all_files)) &
        grepl("DensityU", basename(all_files))
    ]
  }

  if (length(matched) == 0) {
    return(NA_character_)
  }

  return(matched[1])
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
# ES-ALPHA FUNCTIONS
#============================================================

find_beta_file <- function(chr, tauuse) {

  candidate_files <- c(
    paste0(
      output_sub_folder, penalty,
      "chr", chr, "_3traits",
      "warmStart", warmStart,
      "Zscale", Zscale,
      "tauuse", tauuse,
      "singleStart", singleStart,
      ".RData"
    ),
    paste0(
      output_sub_folder, penalty,
      "chr", chr, "_3pop",
      "warmStart", warmStart,
      "Zscale", Zscale,
      "tauuse", tauuse,
      "singleStart", singleStart,
      ".RData"
    ),
    paste0(
      output_sub_folder, penalty,
      "chr", chr, "_3pops",
      "warmStart", warmStart,
      "Zscale", Zscale,
      "tauuse", tauuse,
      "singleStart", singleStart,
      ".RData"
    )
  )

  existing <- candidate_files[file.exists(candidate_files)]

  if (length(existing) > 0) {
    return(existing[1])
  }

  all_files <- list.files(
    output_sub_folder,
    pattern = "\\.RData$",
    full.names = TRUE
  )

  matched <- all_files[
    grepl(paste0("chr", chr), basename(all_files)) &
      grepl(as.character(tauuse), basename(all_files), fixed = TRUE)
  ]

  if (length(matched) == 0) {
    matched <- all_files[
      grepl(paste0("chr", chr), basename(all_files))
    ]
  }

  if (length(matched) == 0) {
    return(NA_character_)
  }

  return(matched[1])
}

compute_es_alpha_from_beta_allchr <- function(
  selected_roww,
  tauuse
) {

  BetaMatrix1_all <- NULL
  BetaMatrix2_all <- NULL
  BetaMatrix3_all <- NULL

  valid_chr_beta <- 0

  for (chr in 1:22) {

    cat("Reading beta matrices for chr", chr, "...\n")

    savefile1 <- find_beta_file(chr = chr, tauuse = tauuse)

    if (is.na(savefile1) || !file.exists(savefile1)) {
      warning("Missing beta file for chr=", chr)
      next
    }

    BetaMatrix1 <- Gen_One_BetaMatrix(savefile1, K, 1)[[1]]
    BetaMatrix2 <- Gen_One_BetaMatrix(savefile1, K, 2)[[1]]
    BetaMatrix3 <- Gen_One_BetaMatrix(savefile1, K, 3)[[1]]

    if (selected_roww > nrow(BetaMatrix1)) {
      stop("selected_roww exceeds number of beta rows for chr ", chr)
    }

    if (is.null(BetaMatrix1_all)) {
      BetaMatrix1_all <- BetaMatrix1
      BetaMatrix2_all <- BetaMatrix2
      BetaMatrix3_all <- BetaMatrix3
    } else {
      BetaMatrix1_all <- cbind(BetaMatrix1_all, BetaMatrix1)
      BetaMatrix2_all <- cbind(BetaMatrix2_all, BetaMatrix2)
      BetaMatrix3_all <- cbind(BetaMatrix3_all, BetaMatrix3)
    }

    valid_chr_beta <- valid_chr_beta + 1
  }

  if (is.null(BetaMatrix1_all)) {
    stop("No valid beta matrices found across chromosomes.")
  }

  cat("\nFinal combined beta matrix dimensions:\n")
  cat(popvec[1], ":", nrow(BetaMatrix1_all), "x", ncol(BetaMatrix1_all), "\n")
  cat(popvec[2], ":", nrow(BetaMatrix2_all), "x", ncol(BetaMatrix2_all), "\n")
  cat(popvec[3], ":", nrow(BetaMatrix3_all), "x", ncol(BetaMatrix3_all), "\n")

  n_alpha_rows <- nrow(BetaMatrix1_all)

  es_alphaMatrix <- matrix(
    0,
    nrow = n_alpha_rows,
    ncol = length(pi_names)
  )

  colnames(es_alphaMatrix) <- pi_names

  for (i in seq_len(n_alpha_rows)) {

    vec1 <- BetaMatrix1_all[i, ]
    vec2 <- BetaMatrix2_all[i, ]
    vec3 <- BetaMatrix3_all[i, ]

    # Correct population order:
    # vec1 = AFR, vec2 = EAS, vec3 = EUR.
    pi_0 <- mean(vec1 == 0 & vec2 == 0 & vec3 == 0)

    pi_1 <- mean(vec1 != 0 & vec2 == 0 & vec3 == 0)
    pi_2 <- mean(vec1 == 0 & vec2 != 0 & vec3 == 0)
    pi_3 <- mean(vec1 == 0 & vec2 == 0 & vec3 != 0)

    pi_12 <- mean(vec1 != 0 & vec2 != 0 & vec3 == 0)
    pi_13 <- mean(vec1 != 0 & vec2 == 0 & vec3 != 0)
    pi_23 <- mean(vec1 == 0 & vec2 != 0 & vec3 != 0)

    pi_123 <- mean(vec1 != 0 & vec2 != 0 & vec3 != 0)

    es_alphaMatrix[i, ] <- c(
      pi_0,
      pi_1,
      pi_2,
      pi_3,
      pi_12,
      pi_13,
      pi_23,
      pi_123
    )
  }

  row_sums <- rowSums(es_alphaMatrix)

  if (max(abs(row_sums - 1), na.rm = TRUE) > 1e-8) {
    warning("Some ES-alpha rows do not sum to 1 before projection.")
  }

  selected_row_es <- es_alphaMatrix[selected_roww, ]
  selected_row_es <- project_simplex(selected_row_es)
  names(selected_row_es) <- pi_names

  list(
    es_alphaMatrix = es_alphaMatrix,
    selected_row = selected_row_es,
    valid_chr_beta = valid_chr_beta,
    total_snps_beta = ncol(BetaMatrix1_all)
  )
}

#============================================================
# SIMPLE CORRECTION FUNCTION
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
  # pi_1 = AFR, pi_2 = EAS, pi_3 = EUR.
  # Shrink each ancestry-specific component proportionally.
  # This avoids artificially forcing the smallest ancestry-specific component,
  # especially EAS/pi_2, to zero.
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
# MAIN REAL-DATA COMPUTATION
#============================================================

rho_idx <- get_rho_index(selected_rho, rho_vec)

cat("\n============================================================\n")
cat("STEP 1: COMPUTE REAL-DATA ES-ALPHA ACROSS ALL CHROMOSOMES\n")
cat("============================================================\n")

es_out <- compute_es_alpha_from_beta_allchr(
  selected_roww = selected_roww,
  tauuse = tauuse
)

es_alphaMatrix <- es_out$es_alphaMatrix
pi_es_alpha <- es_out$selected_row
names(pi_es_alpha) <- pi_names

cat("\nSelected ES-alpha row:\n")
print(pi_es_alpha)

cat("\n============================================================\n")
cat("STEP 2: COMPUTE POSTERIOR SOFT COUNTS FROM REAL DENSITYU FILES\n")
cat("============================================================\n")

pooled_counts <- rep(0, 8)
names(pooled_counts) <- pi_names

total_snps_used <- 0
valid_chr <- 0
density_files_used <- character(0)

for (chr in 1:22) {

  cat("Reading DensityU for chr", chr, "...\n")

  f <- find_density_file(
    chr = chr,
    tauuse = tauuse,
    density_dir = output_sub_folder_uDensity
  )

  if (is.na(f) || !file.exists(f)) {
    warning("Missing DensityU file for chr=", chr)
    next
  }

  if (file.info(f)$size == 0) {
    warning("Empty DensityU file for chr=", chr, ": ", f)
    next
  }

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
    warning("Allrho_DensityU_list missing in: ", f)
    rm(e)
    next
  }

  dens_list <- get("Allrho_DensityU_list", envir = e)

  if (length(dens_list) < rho_idx) {
    warning("Density list shorter than rho_idx for chr=", chr)
    rm(e)
    next
  }

  dens_mat <- dens_list[[rho_idx]]

  if (is.null(dens_mat)) {
    warning("dens_mat is NULL for chr=", chr)
    rm(e)
    next
  }

  if (selected_roww < 1 || selected_roww > nrow(dens_mat)) {
    stop(
      "selected_row out of bounds for chr ",
      chr,
      ". selected_roww=", selected_roww,
      ", nrow(dens_mat)=",
      nrow(dens_mat)
    )
  }

  one <- density_row_to_soft_counts(
    dens_mat[selected_roww, ]
  )

  pooled_counts <- pooled_counts + one$soft_counts
  total_snps_used <- total_snps_used + one$P
  valid_chr <- valid_chr + 1
  density_files_used <- c(density_files_used, f)

  rm(e)
  gc()
}

if (sum(pooled_counts) <= 0) {
  stop("No valid posterior counts were obtained from DensityU files.")
}

pi_final <- project_simplex(pooled_counts)
names(pi_final) <- pi_names

cat("\nRaw posterior soft-count estimate Pi_Final:\n")
print(pi_final)

#============================================================
# BASIC NELDER-MEAD OPTIMIZATION FOR PI
#============================================================

cat("\n============================================================\n")
cat("STEP 3: BASIC NELDER-MEAD CHECK\n")
cat("============================================================\n")

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
  nm_function_count <- unname(as.numeric(nm_from_es$counts["function"]))

} else {

  pi_optimized_nm <- nm_from_final$pi
  nm_start_used <- "pi_final"
  nm_value <- nm_from_final$value
  nm_convergence <- nm_from_final$convergence
  nm_function_count <- unname(as.numeric(nm_from_final$counts["function"]))
}

names(pi_optimized_nm) <- pi_names

cat("Nelder-Mead start used:", nm_start_used, "\n")
cat("Nelder-Mead negative log-likelihood:", nm_value, "\n")
cat("Nelder-Mead convergence:", nm_convergence, "\n")

#============================================================
# DIAGNOSTIC GRIDS WITHOUT TRUE VALUES
#============================================================

cat("\n============================================================\n")
cat("STEP 4: DIAGNOSTIC GRIDS WITHOUT TRUE VALUES\n")
cat("============================================================\n")

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

#============================================================
# FINAL ESTIMATE WITHOUT USING TRUE VALUES
#============================================================

cat("\n============================================================\n")
cat("STEP 5: FINAL REAL-DATA-SAFE ESTIMATE\n")
cat("============================================================\n")

# Secondary sensitivity estimate:
# Simple corrected pi123 rule.
pi_corrected_123 <- correct_pi123_v2(
  pi_vec = pi_final,
  lambda_single = lambda_single_calibrated,
  lambda_zero = lambda_zero_calibrated
)

names(pi_corrected_123) <- pi_names

# Primary final estimate:
# Directional simulation-calibrated method.
directional_calibrated_out <- optimize_pi_nelder_mead_directional(
  soft_counts = pooled_counts,
  start_pi = pi_final,
  lambda_single = lambda_single_directional_calibrated,
  lambda_shared = lambda_shared_directional_calibrated
)

pi_final_optimized <- directional_calibrated_out$pi
pi_final_optimized <- project_simplex(pi_final_optimized)
names(pi_final_optimized) <- pi_names

# Keep this variable for output compatibility.
pi_optimized_directional_nm <- pi_final_optimized
names(pi_optimized_directional_nm) <- pi_names

cat("\nFinal selected estimate Pi_Final_Optimized_NoTrue:\n")
print(pi_final_optimized)

cat("\nSecondary simple corrected pi123 estimate Pi_Best_Corrected_123_NoTrue:\n")
print(pi_corrected_123)

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
    function_count = nm_function_count,
    final_selected = FALSE,
    note = nm_start_used
  ),
  make_candidate_row_no_true(
    method = "simple_corrected_pi123_v2_sensitivity_no_true",
    pi_vec = pi_corrected_123,
    reference_pi = pi_final,
    lambda_single = lambda_single_calibrated,
    lambda_zero = lambda_zero_calibrated,
    objective_value = sum((pi_corrected_123 - pi_final)^2, na.rm = TRUE),
    convergence = 0,
    final_selected = FALSE,
    note = "Secondary sensitivity estimate; lambdas fixed from simulation calibration."
  ),
  make_candidate_row_no_true(
    method = "final_directional_calibrated_no_true",
    pi_vec = pi_final_optimized,
    reference_pi = pi_final,
    lambda_single = lambda_single_directional_calibrated,
    lambda_shared = lambda_shared_directional_calibrated,
    objective_value = directional_calibrated_out$value,
    convergence = directional_calibrated_out$convergence,
    function_count = unname(as.numeric(directional_calibrated_out$counts["function"])),
    final_selected = TRUE,
    note = "Final estimate; directional lambdas fixed from simulation calibration, not selected by true_pi."
  ),
  make_candidate_row_no_true(
    method = "pi_es_alpha",
    pi_vec = pi_es_alpha,
    reference_pi = pi_final,
    objective_value = NA_real_,
    convergence = 0,
    final_selected = FALSE,
    note = "External ES-alpha comparison computed from beta zero/nonzero patterns."
  )
)

rownames(candidate_grid_df) <- NULL

#============================================================
# RESULT TABLE
#============================================================

real_result <- data.frame(
  selected_row = selected_roww,
  selected_rho = selected_rho,
  tauuse = tauuse,
  Parameter = pi_names,
  True_Comparison_Only = rep(NA_real_, length(pi_names)),

  # Raw posterior soft-count estimate.
  Pi_Final = unname(as.numeric(pi_final)),

  # Basic likelihood check.
  Pi_Optimized_NelderMead = unname(as.numeric(pi_optimized_nm)),

  # Directional calibrated estimate.
  # This is the same as Pi_Final_Optimized_NoTrue because directional is final.
  Pi_Directional_Calibrated = unname(as.numeric(pi_optimized_directional_nm)),

  # Simple corrected estimate.
  # Kept as secondary sensitivity only.
  Pi_Best_Corrected_123_NoTrue = unname(as.numeric(pi_corrected_123)),

  # Final estimator.
  Pi_Final_Optimized_NoTrue = unname(as.numeric(pi_final_optimized)),

  # ES-alpha comparison.
  Pi_es_alpha = unname(as.numeric(pi_es_alpha)),

  Final_Estimator_Note =
    "Pi_Final_Optimized_NoTrue is the final directional simulation-calibrated estimate using fixed lambdas, not true_pi.",
  stringsAsFactors = FALSE
)

real_metrics <- data.frame(
  selected_row = selected_roww,
  selected_rho = selected_rho,
  tauuse = tauuse,

  total_snps_used_density = total_snps_used,
  valid_chr_density = valid_chr,
  total_snps_used_beta = es_out$total_snps_beta,
  valid_chr_beta = es_out$valid_chr_beta,

  NelderMead_Start = nm_start_used,
  NelderMead_NegLogLik = nm_value,
  NelderMead_Convergence = nm_convergence,
  NelderMead_Function_Count = nm_function_count,

  Final_Optimized_Method = "final_directional_calibrated_no_true",
  Final_Optimized_Lambda_Single = lambda_single_directional_calibrated,
  Final_Optimized_Lambda_Shared = lambda_shared_directional_calibrated,
  Final_Optimized_Lambda_Zero = NA_real_,
  Final_Optimized_Gamma = NA_real_,
  Final_Optimized_Objective = directional_calibrated_out$value,
  Final_Optimized_Convergence = directional_calibrated_out$convergence,
  Final_Optimized_Function_Count =
    unname(as.numeric(directional_calibrated_out$counts["function"])),

  Final_Parameter_Choice_Note =
    "Final directional lambdas fixed from aggregate simulation calibration; true_pi not used.",

  Simple_Corrected_Lambda_Single = lambda_single_calibrated,
  Simple_Corrected_Lambda_Zero = lambda_zero_calibrated,
  Simple_Corrected_Objective =
    sum((pi_corrected_123 - pi_final)^2, na.rm = TRUE),

  Total_Abs_Change_Final_vs_Raw =
    sum(abs(pi_final_optimized - pi_final), na.rm = TRUE),
  Max_Abs_Change_Final_vs_Raw =
    max(abs(pi_final_optimized - pi_final), na.rm = TRUE),

  Total_Abs_Change_Simple_vs_Raw =
    sum(abs(pi_corrected_123 - pi_final), na.rm = TRUE),
  Max_Abs_Change_Simple_vs_Raw =
    max(abs(pi_corrected_123 - pi_final), na.rm = TRUE),

  Pi123_Raw = unname(as.numeric(pi_final["pi_123"])),
  Pi123_Simple_Corrected_NoTrue = unname(as.numeric(pi_corrected_123["pi_123"])),
  Pi123_Final_Optimized_NoTrue = unname(as.numeric(pi_final_optimized["pi_123"])),

  stringsAsFactors = FALSE
)

rownames(real_result) <- NULL
rownames(real_metrics) <- NULL

#============================================================
# PRINT RESULTS
#============================================================

cat("\n==== FINAL PI COMPARISON TABLE: REAL DATA ====\n")
print(real_result)

cat("\n==== REAL-DATA METRICS / DIAGNOSTICS ====\n")
print(real_metrics)

cat("\n==== FINAL SELECTED CANDIDATE ====\n")
print(
  candidate_grid_df[
    candidate_grid_df$final_selected == TRUE,
    c(
      "method",
      "lambda_single",
      "lambda_shared",
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

cat("\n==== SECONDARY SIMPLE CORRECTION CANDIDATE ====\n")
print(
  candidate_grid_df[
    candidate_grid_df$method == "simple_corrected_pi123_v2_sensitivity_no_true",
    c(
      "method",
      "lambda_single",
      "lambda_shared",
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

#============================================================
# SAVE FINAL OUTPUTS
#============================================================

out_results <- file.path(
  dirOutput,
  "RealData_Final_pi_notrue_results.csv"
)

out_metrics <- file.path(
  dirOutput,
  "RealData_Final_pi_notrue_metrics.csv"
)

out_candidate_summary <- file.path(
  dirOutput,
  "RealData_Final_pi_candidate_summary_NO_TRUE_USED_FOR_SELECTION.csv"
)

out_correction_grid <- file.path(
  dirOutput,
  "RealData_Final_pi_correction_grid_diagnostic_NO_TRUE.csv"
)

out_directional_grid <- file.path(
  dirOutput,
  "RealData_Final_pi_directional_grid_diagnostic_NO_TRUE.csv"
)

out_unsupervised_score_grid <- file.path(
  dirOutput,
  "RealData_Final_pi_unsupervised_lambda_score_grid_NO_TRUE.csv"
)

out_es_alpha_matrix <- file.path(
  dirOutput,
  "RealData_es_alphaMatrix_combined_allchr_corrected.csv"
)

out_es_alpha_selected <- file.path(
  dirOutput,
  paste0(
    "RealData_es_alpha_selected_row",
    selected_roww,
    "_rho",
    selected_rho,
    "_tau",
    tauuse,
    ".csv"
  )
)

out_density_files_used <- file.path(
  dirOutput,
  "RealData_density_files_used.txt"
)

write.csv(real_result, out_results, row.names = FALSE)
write.csv(real_metrics, out_metrics, row.names = FALSE)
write.csv(candidate_grid_df, out_candidate_summary, row.names = FALSE)
write.csv(correction_grid_df, out_correction_grid, row.names = FALSE)
write.csv(directional_grid_df, out_directional_grid, row.names = FALSE)
write.csv(unsupervised_score_grid_df, out_unsupervised_score_grid, row.names = FALSE)
write.csv(es_alphaMatrix, out_es_alpha_matrix, row.names = FALSE)

write.csv(
  as.data.frame(t(pi_es_alpha)),
  out_es_alpha_selected,
  row.names = FALSE
)

writeLines(density_files_used, out_density_files_used)

cat("\nSaved files:\n")
cat(out_results, "\n")
cat(out_metrics, "\n")
cat(out_candidate_summary, "\n")
cat(out_correction_grid, "\n")
cat(out_directional_grid, "\n")
cat(out_unsupervised_score_grid, "\n")
cat(out_es_alpha_matrix, "\n")
cat(out_es_alpha_selected, "\n")
cat(out_density_files_used, "\n")

cat("\nIMPORTANT:\n")
cat("Pi_Final_Optimized_NoTrue is the FINAL directional calibrated estimator.\n")
cat("Pi_Best_Corrected_123_NoTrue is only a secondary simple-correction sensitivity estimate.\n")
cat("true_pi is not used anywhere in this real-data script.\n")
cat("No RMSE/MAE/Bias value is used to choose the final estimate in real data.\n")
cat("Final lambdas were fixed before real-data application based on simulation calibration.\n")