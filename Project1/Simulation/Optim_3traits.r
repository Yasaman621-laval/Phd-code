rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
})

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

dirBase <- "/lustre10/scratch/yatah3/yatah3/simulation/SimuGenotype/three_traits/"

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

density_dir <- paste0(sc_dir, "/Trans2-new-3traits/")

summary_file <- list.files(
  path = sc_dir,
  pattern = paste0(
    "^MaxLogSum_3traits_sim(10|[1-9])_scenario",
    scenarioIndex,
    "\\.csv$"
  ),
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

lambda_grid <- seq(0, 0.8, by = 0.05)

directional_lambda_grid <- expand.grid(
  lambda_single = c(0, 10, 25, 50, 100, 250, 500, 1000),
  lambda_shared = c(0, 10, 25, 50, 100, 250, 500, 1000)
)

lambda_single_final <- 10
lambda_shared_final <- 0

pi_table <- matrix(
  c(
    0.900829767, 0.003272394, 0.003272394, 0.003151734,
    0.012709033, 0.012728922, 0.012679862, 0.051355895,

    0.900978271, 0.003420898, 0.003264438, 0.003410290,
    0.012370921, 0.012230372, 0.012272802, 0.052052008,

    0.901317709, 0.003614483, 0.003630395, 0.003623765,
    0.011429511, 0.011584645, 0.011295592, 0.053503901,

    0.902298896, 0.004061322, 0.004238996, 0.004128944,
    0.009573210, 0.009426032, 0.009245705, 0.057026895,

    0.904047797, 0.004420648, 0.004512138, 0.004333137,
    0.007160018, 0.007064552, 0.007100352, 0.061361358
  ),
  nrow = 5,
  byrow = TRUE
)

colnames(pi_table) <- pi_names

true_pi <- pi_table[scenarioIndex, ]
names(true_pi) <- pi_names

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

  # Penalize single-trait components
  penalty_single <-
    pi_vec["pi_1"] +
    pi_vec["pi_2"] +
    pi_vec["pi_3"]

  # Reward larger pi_123 when lambda_shared > 0
  reward_shared <- -log(pi_vec["pi_123"] + eps)

  value <-
    neg_ll +
    lambda_single * penalty_single +
    lambda_shared * reward_shared

  if (!is.finite(value)) {
    value <- .Machine$double.xmax
  }

  return(value)
}

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

run_directional_lambda_grid <- function(
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

    grid_rows[[length(grid_rows) + 1]] <- data.frame(
      lambda_single = lambda_single,
      lambda_shared = lambda_shared,
      objective_value = opt$value,
      convergence = opt$convergence,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, grid_rows)
}

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
      "_3traits",
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

calc_metrics <- function(est, true) {

  data.frame(
    MAE = mean(abs(est - true), na.rm = TRUE),
    RMSE = sqrt(mean((est - true)^2, na.rm = TRUE))
  )
}

correct_pi123 <- function(pi_vec, lambda = 0.30) {

  pi_new <- pi_vec

  single_names <- c("pi_1", "pi_2", "pi_3")

  single_mass <- sum(pi_vec[single_names], na.rm = TRUE)

  if (!is.finite(single_mass) || single_mass <= 0) {
    return(project_simplex(pi_new))
  }

  move_mass <- lambda * single_mass

  pi_new[single_names] <- pi_vec[single_names] -
    move_mass / length(single_names)

  pi_new["pi_123"] <- pi_vec["pi_123"] + move_mass

  pi_new <- project_simplex(pi_new)
  names(pi_new) <- names(pi_vec)

  return(pi_new)
}

grid_search_lambda <- function(pi_vec, lambda_grid, gamma = 10) {

  results <- list()

  for (lambda in lambda_grid) {

    pi_corr <- correct_pi123(pi_vec, lambda)

    distortion <- sum((pi_corr - pi_vec)^2, na.rm = TRUE)
    gain_123 <- pi_corr["pi_123"]
    score <- gain_123 - gamma * distortion

    results[[length(results) + 1]] <- data.frame(
      lambda = lambda,
      pi_123 = gain_123,
      distortion = distortion,
      score = score
    )
  }

  results_df <- do.call(rbind, results)

  best_row <- results_df[which.max(results_df$score), ]

  list(
    results = results_df,
    best_lambda = best_row$lambda
  )
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
lambda_grid_results_all <- list()
directional_grid_results_all <- list()

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

    load(es_alpha_file)

    if (exists("selected_row") && length(selected_row) == 8) {
      pi_es_alpha <- project_simplex(selected_row)
    } else {
      warning("selected_row missing or invalid inside: ", es_alpha_file)
      pi_es_alpha <- rep(NA_real_, 8)
    }

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
  # DIRECTIONAL NELDER-MEAD OPTIMIZATION FOR PI
  # FINAL CHOICE DOES NOT USE TRUE_PI
  #============================================================

  directional_out <- optimize_pi_nelder_mead_directional(
    soft_counts = pooled_counts,
    start_pi = pi_final,
    lambda_single = lambda_single_final,
    lambda_shared = lambda_shared_final
  )

  pi_optimized_directional_nm <- directional_out$pi
  names(pi_optimized_directional_nm) <- pi_names

  if (directional_out$convergence != 0) {

    warning(
      "Directional Nelder-Mead did not converge for sim ",
      sim,
      ". Falling back to pi_final."
    )

    pi_optimized_directional_nm <- pi_final
  }

  best_lambda_single <- lambda_single_final
  best_lambda_shared <- lambda_shared_final

  directional_nm_value <- directional_out$value
  directional_nm_convergence <- directional_out$convergence
  directional_nm_function_count <- directional_out$counts["function"]

  directional_grid_df <- run_directional_lambda_grid(
    soft_counts = pooled_counts,
    start_pi = pi_final,
    directional_lambda_grid = directional_lambda_grid
  )

  directional_grid_df$sim <- sim
  directional_grid_df$scenarioIndex <- scenarioIndex
  directional_grid_df$scenarioName <- sc_name

  directional_grid_df$Final_Lambda_Used <-
    directional_grid_df$lambda_single == lambda_single_final &
    directional_grid_df$lambda_shared == lambda_shared_final

  directional_grid_results_all[[length(directional_grid_results_all) + 1]] <-
    directional_grid_df

  #============================================================
  # OLD LAMBDA GRID SEARCH KEPT ONLY AS COMPARISON
  #============================================================

  lambda_out <- grid_search_lambda(pi_final, lambda_grid)

  best_lambda <- lambda_out$best_lambda

  pi_corrected_123 <- correct_pi123(
    pi_final,
    lambda = best_lambda
  )

  lambda_grid_df <- lambda_out$results
  lambda_grid_df$sim <- sim
  lambda_grid_df$scenarioIndex <- scenarioIndex
  lambda_grid_df$scenarioName <- sc_name

  lambda_grid_results_all[[length(lambda_grid_results_all) + 1]] <-
    lambda_grid_df

  #============================================================
  # METRICS ONLY FOR EVALUATION
  #============================================================

  metric_final <- calc_metrics(pi_final, true_pi)
  metric_es <- calc_metrics(pi_es_alpha, true_pi)
  metric_corrected_123 <- calc_metrics(pi_corrected_123, true_pi)
  metric_optimized_nm <- calc_metrics(pi_optimized_nm, true_pi)
  metric_optimized_directional_nm <- calc_metrics(
    pi_optimized_directional_nm,
    true_pi
  )

  sim_result <- data.frame(
    sim = sim,
    scenarioIndex = scenarioIndex,
    scenarioName = sc_name,
    selected_row = selected_roww,
    selected_rho = selected_rho,
    tauuse = tauuse,
    Largest_Sum_Of_Logs = largest_logsum,
    Parameter = pi_names,
    True = as.numeric(true_pi),
    Pi_Final = as.numeric(pi_final),
    Pi_Optimized_NelderMead = as.numeric(pi_optimized_nm),
    Pi_Optimized_Directional_NelderMead =
      as.numeric(pi_optimized_directional_nm),
    Pi_Corrected_123 = as.numeric(pi_corrected_123),
    Pi_es_alpha = as.numeric(pi_es_alpha),
    Bias_Final = as.numeric(pi_final - true_pi),
    Bias_Optimized_NelderMead = as.numeric(pi_optimized_nm - true_pi),
    Bias_Optimized_Directional_NelderMead =
      as.numeric(pi_optimized_directional_nm - true_pi),
    Bias_Corrected_123 = as.numeric(pi_corrected_123 - true_pi),
    Bias_es_alpha = as.numeric(pi_es_alpha - true_pi),
    AbsBias_Final = abs(as.numeric(pi_final - true_pi)),
    AbsBias_Optimized_NelderMead =
      abs(as.numeric(pi_optimized_nm - true_pi)),
    AbsBias_Optimized_Directional_NelderMead =
      abs(as.numeric(pi_optimized_directional_nm - true_pi)),
    AbsBias_Corrected_123 =
      abs(as.numeric(pi_corrected_123 - true_pi)),
    AbsBias_es_alpha =
      abs(as.numeric(pi_es_alpha - true_pi)),
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
    best_lambda = best_lambda,
    NelderMead_Start = nm_start_used,
    NelderMead_NegLogLik = nm_value,
    NelderMead_Convergence = nm_convergence,
    NelderMead_Function_Count = as.numeric(nm_function_count),
    Directional_NelderMead_Best_Lambda_Single = best_lambda_single,
    Directional_NelderMead_Best_Lambda_Shared = best_lambda_shared,
    Directional_NelderMead_Objective = directional_nm_value,
    Directional_NelderMead_Convergence = directional_nm_convergence,
    Directional_NelderMead_Function_Count =
      as.numeric(directional_nm_function_count),
    MAE_Final = metric_final$MAE,
    RMSE_Final = metric_final$RMSE,
    MAE_Optimized_NelderMead = metric_optimized_nm$MAE,
    RMSE_Optimized_NelderMead = metric_optimized_nm$RMSE,
    MAE_Optimized_Directional_NelderMead =
      metric_optimized_directional_nm$MAE,
    RMSE_Optimized_Directional_NelderMead =
      metric_optimized_directional_nm$RMSE,
    MAE_Corrected_123 = metric_corrected_123$MAE,
    RMSE_Corrected_123 = metric_corrected_123$RMSE,
    MAE_es_alpha = metric_es$MAE,
    RMSE_es_alpha = metric_es$RMSE,
    stringsAsFactors = FALSE
  )

  cat("\n==== BASIC NELDER-MEAD OPTIMIZATION FOR SIM", sim, "====\n")
  cat("Start used:", nm_start_used, "\n")
  cat("Negative log-likelihood:", nm_value, "\n")
  cat("Convergence code:", nm_convergence, "\n")

  cat("\n==== DIRECTIONAL NELDER-MEAD OPTIMIZATION FOR SIM", sim, "====\n")
  cat("Final lambda_single used:", best_lambda_single, "\n")
  cat("Final lambda_shared used:", best_lambda_shared, "\n")
  cat("Directional objective value:", directional_nm_value, "\n")
  cat("Convergence code:", directional_nm_convergence, "\n")

  print(directional_grid_df)

  cat("\n==== LAMBDA GRID RESULTS FOR SIM", sim, "====\n")
  print(lambda_grid_df[, c("lambda", "pi_123", "distortion", "score")])

  cat("\nBest lambda for sim", sim, "=", best_lambda, "\n")

  cat("\n==== FINAL PI COMPARISON FOR SIM", sim, "====\n")
  print(
    sim_result[
      ,
      c(
        "Parameter",
        "True",
        "Pi_Final",
        "Pi_Optimized_NelderMead",
        "Pi_Optimized_Directional_NelderMead",
        "Pi_Corrected_123",
        "Pi_es_alpha",
        "Bias_Final",
        "Bias_Optimized_NelderMead",
        "Bias_Optimized_Directional_NelderMead",
        "Bias_Corrected_123",
        "Bias_es_alpha"
      )
    ]
  )

  cat("\n==== ERROR SUMMARY FOR SIM", sim, "====\n")
  print(sim_metrics)

  results_all[[length(results_all) + 1]] <- sim_result
  metrics_all[[length(metrics_all) + 1]] <- sim_metrics
}

#============================================================
# SAVE FINAL OUTPUTS
#============================================================

if (length(results_all) == 0) {
  stop("No simulation produced valid results.")
}

final_results <- do.call(rbind, results_all)
final_metrics <- do.call(rbind, metrics_all)
final_lambda_grid <- do.call(rbind, lambda_grid_results_all)
final_directional_grid <- do.call(rbind, directional_grid_results_all)

out_results <- paste0(
  sc_dir,
  "/Final_pi_posterior_count_directional_neldermead_results.csv"
)

out_metrics <- paste0(
  sc_dir,
  "/Final_pi_posterior_count_directional_neldermead_metrics.csv"
)

out_lambda <- paste0(
  sc_dir,
  "/Final_pi_lambda_grid_results.csv"
)

out_directional_grid <- paste0(
  sc_dir,
  "/Final_pi_directional_neldermead_lambda_grid_results.csv"
)

write.csv(final_results, out_results, row.names = FALSE)
write.csv(final_metrics, out_metrics, row.names = FALSE)
write.csv(final_lambda_grid, out_lambda, row.names = FALSE)
write.csv(final_directional_grid, out_directional_grid, row.names = FALSE)

cat("\nSaved files:\n")
cat(out_results, "\n")
cat(out_metrics, "\n")
cat(out_lambda, "\n")
cat(out_directional_grid, "\n")