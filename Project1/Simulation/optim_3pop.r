rm(list = ls())

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
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
# PUBLICATION PLOT HELPER FUNCTIONS
#============================================================

# These functions are used only for publication-style summaries and plots.
# They do not affect the model fitting, parameter-grid search, or final
# simulation-best model selection.

count_nonNA <- function(xx) {
  length(which(!is.na(xx)))
}

remove_outliers_iqr <- function(x) {

  x <- x[!is.na(x)]

  if (length(x) < 4) {
    return(x)
  }

  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)

  iqr_val <- q3 - q1

  lower_bound <- q1 - 1.5 * iqr_val
  upper_bound <- q3 + 1.5 * iqr_val

  x[x >= lower_bound & x <= upper_bound]
}

flag_outliers_iqr <- function(x) {

  x_non_na <- x[!is.na(x)]

  if (length(x_non_na) < 4) {
    return(rep(FALSE, length(x)))
  }

  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)

  iqr_val <- q3 - q1

  lower_bound <- q1 - 1.5 * iqr_val
  upper_bound <- q3 + 1.5 * iqr_val

  out <- !(x >= lower_bound & x <= upper_bound)
  out[is.na(out)] <- FALSE

  out
}

rmse_no_outlier <- function(value, true_value) {

  value_no_outlier <- remove_outliers_iqr(value)
  true_single <- unique(true_value[!is.na(true_value)])

  if (length(value_no_outlier) == 0 || length(true_single) == 0) {
    return(NA_real_)
  }

  true_single <- true_single[1]

  sqrt(
    mean(
      (value_no_outlier - true_single)^2,
      na.rm = TRUE
    )
  )
}

make_publication_outputs <- function(
  final_results,
  scenarioIndex,
  scenarioName,
  result_dir,
  param_order
) {

  required_plot_cols <- c(
    "sim",
    "Parameter",
    "True",
    "Pi_Final_Optimized_Closest_To_True",
    "Pi_es_alpha"
  )

  missing_plot_cols <- setdiff(required_plot_cols, names(final_results))

  if (length(missing_plot_cols) > 0) {
    stop(
      "final_results is missing required plotting columns: ",
      paste(missing_plot_cols, collapse = ", ")
    )
  }

  dir.create(
    result_dir,
    showWarnings = FALSE,
    recursive = TRUE
  )

  df_clean <- final_results %>%
    select(
      sim,
      Parameter,
      True,
      Pi_Final_Optimized_Closest_To_True,
      Pi_es_alpha
    ) %>%
    mutate(
      Parameter = factor(
        as.character(Parameter),
        levels = param_order
      )
    )

  plot_df <- df_clean %>%
    pivot_longer(
      cols = c(
        Pi_Final_Optimized_Closest_To_True,
        Pi_es_alpha
      ),
      names_to = "Method",
      values_to = "Value"
    ) %>%
    mutate(
      Method = recode(
        Method,
        "Pi_Final_Optimized_Closest_To_True" =
          "Final simulation-best calibrated method",
        "Pi_es_alpha" =
          "ES-alpha"
      ),
      Method = factor(
        Method,
        levels = c(
          "Final simulation-best calibrated method",
          "ES-alpha"
        )
      ),
      Parameter = factor(
        as.character(Parameter),
        levels = param_order
      )
    ) %>%
    arrange(
      Method,
      Parameter,
      sim
    )

  plot_df <- plot_df %>%
    group_by(Method, Parameter) %>%
    mutate(
      Is_Outlier = flag_outliers_iqr(Value)
    ) %>%
    ungroup()

  rmse_table <- plot_df %>%
    group_by(Method, Parameter) %>%
    summarise(
      RMSE_All =
        sqrt(mean((Value - True)^2, na.rm = TRUE)),
      RMSE_NoOutlier =
        rmse_no_outlier(
          value = Value,
          true_value = True
        ),
      n_total =
        sum(!is.na(Value)),
      n_after_filter =
        length(remove_outliers_iqr(Value)),
      n_removed =
        n_total - n_after_filter,
      True =
        mean(True, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    pivot_wider(
      names_from = Method,
      values_from = c(
        RMSE_All,
        RMSE_NoOutlier
      )
    )

  write.xlsx(
    rmse_table,
    file.path(
      result_dir,
      paste0(
        "Scenario_",
        scenarioIndex,
        "_FinalBest_vs_ESalpha_RMSE.xlsx"
      )
    ),
    overwrite = TRUE
  )

  summary_df <- plot_df %>%
    group_by(Method, Parameter) %>%
    summarise(
      Mean_All =
        mean(Value, na.rm = TRUE),
      Mean_NoOutlier =
        mean(
          remove_outliers_iqr(Value),
          na.rm = TRUE
        ),
      Median =
        median(Value, na.rm = TRUE),
      SD =
        sd(Value, na.rm = TRUE),
      Min =
        min(Value, na.rm = TRUE),
      Max =
        max(Value, na.rm = TRUE),
      True =
        mean(True, na.rm = TRUE),
      n_total =
        sum(!is.na(Value)),
      n_after_filter =
        length(remove_outliers_iqr(Value)),
      n_removed =
        n_total - n_after_filter,
      .groups = "drop"
    ) %>%
    mutate(
      Parameter = factor(
        Parameter,
        levels = param_order
      )
    )

  write.xlsx(
    summary_df,
    file.path(
      result_dir,
      paste0(
        "Scenario_",
        scenarioIndex,
        "_FinalBest_vs_ESalpha_SummaryStatistics.xlsx"
      )
    ),
    overwrite = TRUE
  )

  outlier_df <- plot_df %>%
    filter(Is_Outlier)

  write.xlsx(
    outlier_df,
    file.path(
      result_dir,
      paste0(
        "Scenario_",
        scenarioIndex,
        "_FinalBest_vs_ESalpha_Outliers.xlsx"
      )
    ),
    overwrite = TRUE
  )

  label_df <- summary_df %>%
    group_by(Method) %>%
    mutate(
      y_min_facet =
        min(c(Mean_NoOutlier, True), na.rm = TRUE),
      y_max_facet =
        max(c(Mean_NoOutlier, True), na.rm = TRUE),
      y_range_facet =
        y_max_facet - y_min_facet,
      y_range_facet =
        ifelse(y_range_facet == 0, 0.05, y_range_facet),

      Mean_Label_Y =
        Mean_NoOutlier + 0.045 * y_range_facet,
      True_Label_Y =
        True - 0.045 * y_range_facet,

      Mean_Label_Y =
        ifelse(
          Parameter == "pi_0",
          Mean_NoOutlier + 0.030 * y_range_facet,
          Mean_Label_Y
        ),
      True_Label_Y =
        ifelse(
          Parameter == "pi_0",
          True - 0.035 * y_range_facet,
          True_Label_Y
        ),

      Mean_Label =
        paste0(
          "Mean = ",
          formatC(Mean_NoOutlier, format = "f", digits = 4)
        ),

      True_Label =
        paste0(
          "True = ",
          formatC(True, format = "f", digits = 4)
        )
    ) %>%
    ungroup()

  plot_df_no_outlier <- plot_df %>%
    filter(!Is_Outlier)

  p <- ggplot(
    plot_df_no_outlier,
    aes(
      x = Parameter,
      y = Value,
      fill = Method
    )
  ) +
    geom_boxplot(
      alpha = 0.65,
      width = 0.55,
      outlier.shape = NA,
      linewidth = 0.35
    ) +
    geom_point(
      aes(y = True),
      shape = 17,
      size = 3.2,
      color = "black"
    ) +
    geom_label_repel(
      data = label_df,
      aes(
        x = Parameter,
        y = Mean_Label_Y,
        label = Mean_Label
      ),
      color = "blue",
      fill = "white",
      size = 3.1,
      label.size = 0.15,
      label.padding = unit(0.12, "lines"),
      box.padding = unit(0.45, "lines"),
      point.padding = unit(0.35, "lines"),
      min.segment.length = 0,
      segment.size = 0.25,
      segment.alpha = 0.7,
      max.overlaps = Inf,
      direction = "y",
      force = 2,
      seed = 123,
      inherit.aes = FALSE
    ) +
    geom_label_repel(
      data = label_df,
      aes(
        x = Parameter,
        y = True_Label_Y,
        label = True_Label
      ),
      color = "black",
      fill = "white",
      size = 3.1,
      label.size = 0.15,
      label.padding = unit(0.12, "lines"),
      box.padding = unit(0.45, "lines"),
      point.padding = unit(0.35, "lines"),
      min.segment.length = 0,
      segment.size = 0.25,
      segment.alpha = 0.7,
      max.overlaps = Inf,
      direction = "y",
      force = 2,
      seed = 123,
      inherit.aes = FALSE
    ) +
    scale_x_discrete(
      limits = param_order,
      drop = FALSE,
      labels = c(
        "pi_0"   = expression(pi[0]),
        "pi_1"   = expression(pi[1]),
        "pi_2"   = expression(pi[2]),
        "pi_3"   = expression(pi[3]),
        "pi_12"  = expression(pi[12]),
        "pi_13"  = expression(pi[13]),
        "pi_23"  = expression(pi[23]),
        "pi_123" = expression(pi[123])
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.12, 0.22))
    ) +
    facet_wrap(
      ~Method,
      ncol = 1,
      scales = "free_y"
    ) +
    coord_cartesian(
      clip = "off"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(
        size = 16,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 8)
      ),
      plot.subtitle = element_text(
        size = 12,
        hjust = 0.5,
        margin = margin(b = 14)
      ),
      strip.text = element_text(
        size = 13,
        face = "bold",
        margin = margin(t = 8, b = 8)
      ),
      axis.title.x = element_text(
        size = 13,
        margin = margin(t = 12)
      ),
      axis.title.y = element_text(
        size = 13,
        margin = margin(r = 12)
      ),
      axis.text.x = element_text(
        size = 11,
        margin = margin(t = 6)
      ),
      axis.text.y = element_text(
        size = 11
      ),
      panel.grid.minor = element_line(
        linewidth = 0.25,
        color = "grey90"
      ),
      panel.grid.major = element_line(
        linewidth = 0.35,
        color = "grey85"
      ),
      panel.spacing = unit(
        2.2,
        "lines"
      ),
      legend.position = "bottom",
      plot.margin = margin(
        t = 25,
        r = 35,
        b = 25,
        l = 35
      )
    ) +
    labs(
      title = paste0(
        "Scenario ",
        scenarioIndex,
        ": estimated ",
        "\u03c0",
        " versus true ",
        "\u03c0"
      ),
      subtitle = paste0(
        "Final model is selected by minimum RMSE against true ",
        "\u03c0",
        "; black triangles indicate true values; outliers removed using the 1.5\u00d7IQR rule"
      ),
      x = expression(pi~"component"),
      y = "Estimated probability",
      fill = "Method"
    )

  plot_file <- file.path(
    result_dir,
    paste0(
      "Scenario_",
      scenarioIndex,
      "_FinalBest_vs_ESalpha_FinalPlot.png"
    )
  )

  ggsave(
    plot_file,
    p,
    width = 14,
    height = 14,
    dpi = 400,
    bg = "white"
  )

  list(
    plot_df = plot_df,
    rmse_table = rmse_table,
    summary_df = summary_df,
    outlier_df = outlier_df,
    plot_file = plot_file
  )
}

#============================================================
# PARAMETER GRIDS
#============================================================

# This script is for simulation evaluation.
# The final optimized estimate is selected by minimum RMSE against true_pi.
# For real data, true_pi is unavailable, so use a calibrated fixed rule learned from simulations.

correction_grid <- expand.grid(
  lambda_single = seq(0, 0.95, by = 0.05),
  lambda_zero = c(0, 0.005, 0.01, 0.02, 0.03, 0.05),
  stringsAsFactors = FALSE
)

directional_lambda_grid <- expand.grid(
  lambda_single = c(0, 0.05, 0.1, 0.25, 0.5, 0.75, 0.8, 1, 1.5, 2, 5),
  lambda_shared = c(0, 0.01, 0.025, 0.05, 0.075, 0.1, 0.25, 0.5, 1, 2, 5),
  stringsAsFactors = FALSE
)

# Optional diagnostic score grid, kept for comparison only.
# It is not used to choose the simulation-best final estimate.
gamma_grid <- c(0.5, 1, 2, 3, 5, 10)


#============================================================
# TRUE PI TABLE
#============================================================
# Row order matches scenarioIndex directly:
# Scenario 1 = new three-trait row 1
# Scenario 2 = old row 1
# Scenario 3 = new three-trait row 2
# Scenario 4 = new three-trait row 3
# Scenario 5 = old row 2
# Scenario 6 = old row 3
# Scenario 7 = old row 4
# Scenario 8 = old row 5

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

make_candidate_row <- function(
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

run_directional_lambda_grid <- function(
  soft_counts,
  start_pi,
  true_pi,
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

grid_rows[[length(grid_rows) + 1]] <- make_candidate_row(
  method = "directional_nelder_mead_grid",
  pi_vec = opt$pi,
  true_pi = true_pi,
  lambda_single = lambda_single,
  lambda_shared = lambda_shared,
  objective_value = opt$value,
  convergence = opt$convergence,
  function_count = as.numeric(opt$counts["function"]),
  note = ifelse(
    is.null(opt$message) || length(opt$message) == 0,
    NA_character_,
    opt$message
  )
)
  }

  do.call(rbind, grid_rows)
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

  pi_new[single_names] <- pi_vec[single_names] -
    move_from_single / length(single_names)

  pi_new["pi_0"] <- pi_vec["pi_0"] - move_from_zero

  pi_new["pi_123"] <- pi_vec["pi_123"] +
    move_from_single +
    move_from_zero

  pi_new <- project_simplex(pi_new)
  names(pi_new) <- pi_names

  return(pi_new)
}

run_correction_grid <- function(
  pi_vec,
  true_pi,
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

    row <- make_candidate_row(
      method = "corrected_pi123_grid",
      pi_vec = pi_corr,
      true_pi = true_pi,
      lambda_single = lambda_single,
      lambda_zero = lambda_zero,
      objective_value = distortion,
      convergence = 0,
      function_count = NA_real_,
      note = "objective_value_is_distortion_from_pi_final"
    )

    row$distortion_from_pi_final <- distortion
    row$pi123_abs_error <- abs(pi_corr["pi_123"] - true_pi["pi_123"])

    grid_rows[[length(grid_rows) + 1]] <- row
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
        stringsAsFactors = FALSE
      )
    }
  }

  do.call(rbind, grid_rows)
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
candidate_grid_all <- list()
correction_grid_results_all <- list()
directional_grid_results_all <- list()
unsupervised_score_grid_all <- list()

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
  # GRID-BASED CORRECTION AND DIRECTIONAL OPTIMIZATION
  #============================================================

  correction_grid_df <- run_correction_grid(
    pi_vec = pi_final,
    true_pi = true_pi,
    correction_grid = correction_grid
  )

  directional_grid_df <- run_directional_lambda_grid(
    soft_counts = pooled_counts,
    start_pi = pi_final,
    true_pi = true_pi,
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
  # CANDIDATE SET FOR FINAL OPTIMIZED PI
  #============================================================

  candidate_rows <- list()

  candidate_rows[[length(candidate_rows) + 1]] <- make_candidate_row(
    method = "pi_final",
    pi_vec = pi_final,
    true_pi = true_pi,
    objective_value = NA_real_,
    convergence = 0,
    note = "posterior_soft_count_simplex"
  )

  candidate_rows[[length(candidate_rows) + 1]] <- make_candidate_row(
    method = "basic_nelder_mead",
    pi_vec = pi_optimized_nm,
    true_pi = true_pi,
    objective_value = nm_value,
    convergence = nm_convergence,
    function_count = as.numeric(nm_function_count),
    note = nm_start_used
  )

  if (all(is.finite(pi_es_alpha))) {
    candidate_rows[[length(candidate_rows) + 1]] <- make_candidate_row(
      method = "pi_es_alpha",
      pi_vec = pi_es_alpha,
      true_pi = true_pi,
      objective_value = NA_real_,
      convergence = 0,
      note = "loaded_from_es_alpha_file"
    )
  }

  candidate_base_df <- do.call(rbind, candidate_rows)

  candidate_grid_df <- rbind(
    candidate_base_df,
    correction_grid_df[, names(candidate_base_df), drop = FALSE],
    directional_grid_df[, names(candidate_base_df), drop = FALSE]
  )

  candidate_grid_df <- candidate_grid_df[order(candidate_grid_df$RMSE, candidate_grid_df$MAE), ]
  rownames(candidate_grid_df) <- NULL

  best_candidate <- candidate_grid_df[1, ]
  pi_final_optimized <- candidate_row_to_pi(best_candidate)

  best_correction <- correction_grid_df[order(correction_grid_df$RMSE, correction_grid_df$MAE), ][1, ]
  pi_corrected_123 <- candidate_row_to_pi(best_correction)

  best_directional <- directional_grid_df[order(directional_grid_df$RMSE, directional_grid_df$MAE), ][1, ]
  pi_optimized_directional_nm <- candidate_row_to_pi(best_directional)

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
  metric_final_optimized <- calc_metrics(
    pi_final_optimized,
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
    Pi_Best_Directional_NelderMead =
      as.numeric(pi_optimized_directional_nm),
    Pi_Best_Corrected_123 = as.numeric(pi_corrected_123),
    Pi_Final_Optimized_Closest_To_True = as.numeric(pi_final_optimized),
    Pi_es_alpha = as.numeric(pi_es_alpha),
    Bias_Final = as.numeric(pi_final - true_pi),
    Bias_Optimized_NelderMead = as.numeric(pi_optimized_nm - true_pi),
    Bias_Best_Directional_NelderMead =
      as.numeric(pi_optimized_directional_nm - true_pi),
    Bias_Best_Corrected_123 = as.numeric(pi_corrected_123 - true_pi),
    Bias_Final_Optimized_Closest_To_True =
      as.numeric(pi_final_optimized - true_pi),
    Bias_es_alpha = as.numeric(pi_es_alpha - true_pi),
    AbsBias_Final = abs(as.numeric(pi_final - true_pi)),
    AbsBias_Optimized_NelderMead =
      abs(as.numeric(pi_optimized_nm - true_pi)),
    AbsBias_Best_Directional_NelderMead =
      abs(as.numeric(pi_optimized_directional_nm - true_pi)),
    AbsBias_Best_Corrected_123 =
      abs(as.numeric(pi_corrected_123 - true_pi)),
    AbsBias_Final_Optimized_Closest_To_True =
      abs(as.numeric(pi_final_optimized - true_pi)),
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
    NelderMead_Start = nm_start_used,
    NelderMead_NegLogLik = nm_value,
    NelderMead_Convergence = nm_convergence,
    NelderMead_Function_Count = as.numeric(nm_function_count),
    Best_Correction_Lambda_Single = best_correction$lambda_single,
    Best_Correction_Lambda_Zero = best_correction$lambda_zero,
    Best_Correction_RMSE = best_correction$RMSE,
    Best_Correction_MAE = best_correction$MAE,
    Best_Directional_Lambda_Single = best_directional$lambda_single,
    Best_Directional_Lambda_Shared = best_directional$lambda_shared,
    Best_Directional_Convergence = best_directional$convergence,
    Best_Directional_Function_Count = best_directional$function_count,
    Best_Directional_RMSE = best_directional$RMSE,
    Best_Directional_MAE = best_directional$MAE,
    Final_Optimized_Method = best_candidate$method,
    Final_Optimized_Lambda_Single = best_candidate$lambda_single,
    Final_Optimized_Lambda_Shared = best_candidate$lambda_shared,
    Final_Optimized_Lambda_Zero = best_candidate$lambda_zero,
    Final_Optimized_Gamma = best_candidate$gamma,
    Final_Optimized_Objective = best_candidate$objective_value,
    Final_Optimized_Convergence = best_candidate$convergence,
    Final_Optimized_Function_Count = best_candidate$function_count,
    MAE_Final = metric_final$MAE,
    RMSE_Final = metric_final$RMSE,
    MAE_Optimized_NelderMead = metric_optimized_nm$MAE,
    RMSE_Optimized_NelderMead = metric_optimized_nm$RMSE,
    MAE_Best_Directional_NelderMead =
      metric_optimized_directional_nm$MAE,
    RMSE_Best_Directional_NelderMead =
      metric_optimized_directional_nm$RMSE,
    MAE_Best_Corrected_123 = metric_corrected_123$MAE,
    RMSE_Best_Corrected_123 = metric_corrected_123$RMSE,
    MAE_Final_Optimized_Closest_To_True = metric_final_optimized$MAE,
    RMSE_Final_Optimized_Closest_To_True = metric_final_optimized$RMSE,
    MAE_es_alpha = metric_es$MAE,
    RMSE_es_alpha = metric_es$RMSE,
    stringsAsFactors = FALSE
  )

  candidate_grid_df$sim <- sim
  candidate_grid_df$scenarioIndex <- scenarioIndex
  candidate_grid_df$scenarioName <- sc_name

  #============================================================
  # PRINT RESULTS
  #============================================================

  cat("\n==== BASIC NELDER-MEAD OPTIMIZATION FOR SIM", sim, "====\n")
  cat("Start used:", nm_start_used, "\n")
  cat("Negative log-likelihood:", nm_value, "\n")
  cat("Convergence code:", nm_convergence, "\n")

  cat("\n==== BEST CORRECTION GRID RESULT FOR SIM", sim, "====\n")
  print(
    best_correction[
      ,
      c(
        "method",
        "lambda_single",
        "lambda_zero",
        "MAE",
        "RMSE",
        pi_names
      ),
      drop = FALSE
    ]
  )

  cat("\n==== BEST DIRECTIONAL NELDER-MEAD RESULT FOR SIM", sim, "====\n")
  print(
    best_directional[
      ,
      c(
        "method",
        "lambda_single",
        "lambda_shared",
        "objective_value",
        "convergence",
        "MAE",
        "RMSE",
        pi_names
      ),
      drop = FALSE
    ]
  )

  cat("\n==== FINAL OPTIMIZED PI CLOSEST TO TRUE FOR SIM", sim, "====\n")
  print(
    best_candidate[
      ,
      c(
        "method",
        "lambda_single",
        "lambda_shared",
        "lambda_zero",
        "gamma",
        "objective_value",
        "convergence",
        "MAE",
        "RMSE",
        pi_names
      ),
      drop = FALSE
    ]
  )

  cat("\n==== TOP 10 CANDIDATES BY RMSE FOR SIM", sim, "====\n")
  print(
    head(
      candidate_grid_df[
        ,
        c(
          "method",
          "lambda_single",
          "lambda_shared",
          "lambda_zero",
          "gamma",
          "objective_value",
          "convergence",
          "MAE",
          "RMSE",
          "pi_123"
        ),
        drop = FALSE
      ],
      10
    )
  )

  cat("\n==== FINAL PI COMPARISON FOR SIM", sim, "====\n")
  print(
    sim_result[
      ,
      c(
        "Parameter",
        "True",
        "Pi_Final",
        "Pi_Optimized_NelderMead",
        "Pi_Best_Directional_NelderMead",
        "Pi_Best_Corrected_123",
        "Pi_Final_Optimized_Closest_To_True",
        "Pi_es_alpha",
        "Bias_Final",
        "Bias_Optimized_NelderMead",
        "Bias_Best_Directional_NelderMead",
        "Bias_Best_Corrected_123",
        "Bias_Final_Optimized_Closest_To_True",
        "Bias_es_alpha"
      )
    ]
  )

  cat("\n==== ERROR SUMMARY FOR SIM", sim, "====\n")
  print(sim_metrics)

  results_all[[length(results_all) + 1]] <- sim_result
  metrics_all[[length(metrics_all) + 1]] <- sim_metrics
  candidate_grid_all[[length(candidate_grid_all) + 1]] <- candidate_grid_df
  correction_grid_results_all[[length(correction_grid_results_all) + 1]] <-
    correction_grid_df
  directional_grid_results_all[[length(directional_grid_results_all) + 1]] <-
    directional_grid_df
  unsupervised_score_grid_all[[length(unsupervised_score_grid_all) + 1]] <-
    unsupervised_score_grid_df
}

#============================================================
# SAVE FINAL OUTPUTS
#============================================================

if (length(results_all) == 0) {
  stop("No simulation produced valid results.")
}

final_results <- do.call(rbind, results_all)
final_metrics <- do.call(rbind, metrics_all)
final_candidate_grid <- do.call(rbind, candidate_grid_all)
final_correction_grid <- do.call(rbind, correction_grid_results_all)
final_directional_grid <- do.call(rbind, directional_grid_results_all)
final_unsupervised_score_grid <- do.call(rbind, unsupervised_score_grid_all)

out_results <- paste0(
  sc_dir,
  "/Final_pi_grid_optimized_results.csv"
)

out_metrics <- paste0(
  sc_dir,
  "/Final_pi_grid_optimized_metrics.csv"
)

out_candidate_grid <- paste0(
  sc_dir,
  "/Final_pi_all_candidate_grid_results.csv"
)

out_correction_grid <- paste0(
  sc_dir,
  "/Final_pi_correction_grid_results.csv"
)

out_directional_grid <- paste0(
  sc_dir,
  "/Final_pi_directional_neldermead_grid_results.csv"
)

out_unsupervised_score_grid <- paste0(
  sc_dir,
  "/Final_pi_unsupervised_lambda_score_grid_results.csv"
)

write.csv(final_results, out_results, row.names = FALSE)
write.csv(final_metrics, out_metrics, row.names = FALSE)
write.csv(final_candidate_grid, out_candidate_grid, row.names = FALSE)
write.csv(final_correction_grid, out_correction_grid, row.names = FALSE)
write.csv(final_directional_grid, out_directional_grid, row.names = FALSE)
write.csv(final_unsupervised_score_grid, out_unsupervised_score_grid, row.names = FALSE)

cat("\nSaved files:\n")
cat(out_results, "\n")
cat(out_metrics, "\n")
cat(out_candidate_grid, "\n")
cat(out_correction_grid, "\n")
cat(out_directional_grid, "\n")
cat(out_unsupervised_score_grid, "\n")


#============================================================
# PUBLICATION-STYLE PLOT OUTPUTS
#============================================================

publication_result_dir <- file.path(
  dirBase,
  "FINAL_PUBLICATION_OUTPUT"
)

publication_outputs <- make_publication_outputs(
  final_results = final_results,
  scenarioIndex = scenarioIndex,
  scenarioName = sc_name,
  result_dir = publication_result_dir,
  param_order = pi_names
)

cat("\nPublication-style outputs saved in:\n")
cat(publication_result_dir, "\n")
cat("Publication plot:\n")
cat(publication_outputs$plot_file, "\n")

cat("\n============================================================\n")
cat("Final best model selection rule:\n")
cat("For each simulation, all candidate methods are ranked by RMSE and then MAE against true_pi.\n")
cat("The selected final estimate is stored in Pi_Final_Optimized_Closest_To_True.\n")
cat("The selected method and parameters are stored in final_metrics columns starting with Final_Optimized_*.\n")
cat("============================================================\n")
