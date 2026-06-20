#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(SummaryLasso)
  library(dplyr)
  library(data.table)
  library(stringr)
})

## ============================================================
## 1. Parse job index (SLURM array ID) & read input.txt
## ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("? Missing argument: SLURM_ARRAY_TASK_ID")
jset <- as.numeric(args[1])

file.rjobs <- "/lustre10/scratch/yatah3/yatah3/simulation/project1/step3/"
inputs     <- read.table(paste0(file.rjobs, "input.txt"),
                         header = TRUE, as.is = TRUE)

## Expected columns: scenarioIndex, sim  (same style as before)
scenarioIndex <- inputs[jset, 1]
sim           <- inputs[jset, 2]

cat(">>> Post-processing script (3 traits, 1 pop: EUR)\n")
cat(">>> sim =", sim, " scenarioIndex =", scenarioIndex, "\n\n")

## ============================================================
## 2. Scenario detection (same logic as step 1 script, but three_traits)
## ============================================================
dirBase <- "/lustre10/scratch/yatah3/yatah3/simulation/SimuGenotype/three_traits/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) {
  stop("No scenario directories found under: ", dirBase)
}
if (scenarioIndex > length(scenario_dirs)) {
  stop("scenarioIndex exceeds available scenarios (",
       length(scenario_dirs), ")")
}

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)
cat(">>> Using scenario directory:", sc_dir, "\n\n")

## ============================================================
## 3. General settings (consistent with 3-traits step1)
## ============================================================
savename  <- "Trans2-new-3traits"
penalty   <- "RealmixLOG"

## one-population, three-traits case:
popvec    <- "EUR"

K         <- 3
epsilon   <- 1e-8
numsnp    <- 400
rho_vec   <- c(seq(0, 0.9, 0.1), 0.95)
warmStart <- 1
Zscale    <- 1
singleStart <- 1

## Directory where step 1 saved the u-density files
## (must match step1: Trans2-new-3traits + "_3traits" in filename)
output_sub_folder_uDensity <- file.path(sc_dir, savename)

## Directory to save combined matrices / summary (optional)
dirtemp      <- Sys.getenv("SLURM_TMPDIR")
if (dirtemp == "" || is.na(dirtemp)) {
  ## Fallback: store in scenario dir if SLURM_TMPDIR not set
  dirtemp <- file.path(sc_dir, "tmp_step3")
}
output_sub_folder_Beta <- paste0(dirtemp, "/Allchr_BetaMatrices/")
dir.create(output_sub_folder_Beta, recursive = TRUE, showWarnings = FALSE)

## ============================================================
## 4. Load Tau matrix (AbsTauvec) for this scenario
##    (from step1: Taufile <- paste0(sc_dir,"/",popvec,"3traits_Tau_info_v2.RData"))
## ============================================================
Taufile <- file.path(
  sc_dir,
  paste0(popvec, "3traits_Tau_info_v2.RData")
)

if (!file.exists(Taufile)) {
  stop("Taufile missing: ", Taufile)
}
load(Taufile)  ## must load AbsTauvec
if (!exists("AbsTauvec")) {
  stop("Object 'AbsTauvec' not found in Taufile: ", Taufile)
}

cat(">>> Taufile loaded. dim(AbsTauvec) = ",
    paste(dim(AbsTauvec), collapse = " x "), "\n\n")

## ============================================================
## 5. Containers to track global maximum
## ============================================================
largest_values  <- numeric(0)
largest_indices <- numeric(0)
largest_rhos    <- numeric(0)
largest_tauuse  <- numeric(0)

## ============================================================
## 6. Main loop over Tau columns (tt) and rho
##    Structure kept same as original 3-pop code
## ============================================================
nTauCols <- ncol(AbsTauvec)

for (tt in seq_len(nTauCols)) {

  tauuse <- AbsTauvec[sim, tt]
  cat(">>> Tau column tt =", tt, " tauuse =", tauuse, "\n")

  ## Skip if tauuse is NA (safety)
  if (is.na(tauuse)) {
    cat("    - tauuse is NA, skipping this column.\n")
    next
  }

  ## List to store, for this tauuse, the full AllBeta across all chrs
  ## for each rho index rr
  Allrho_ZMatrix_list <- vector("list", length(rho_vec))

  ## ----------------------------------
  ## 6a. For each rho, combine across chromosomes
  ## ----------------------------------
  for (rr in seq_along(rho_vec)) {
    rho <- rho_vec[rr]
    cat("    > rho index rr =", rr, " rho =", rho, "\n")

    AllBeta <- NULL

    for (chr in 1:22) {

      ## Must match step1 filename pattern:
      ## penalty, "_chr", chr, "_3traits", "_warmStart", warmStart,
      ## "_sim", sim, "_Zscale", Zscale, "_singleStart", singleStart,
      ## "_tauuse", tauuse, "_DensityU.RData"
      savefile1 <- file.path(
        output_sub_folder_uDensity,
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

      if (!file.exists(savefile1) || file.info(savefile1)$size <= 0) {
        cat("       - Missing file for chr =", chr,
            ", tauuse =", tauuse, " (", savefile1, ")\n")
        next
      }

      ## Load Allrho_WBMatrix_list for this (chr, tauuse)
      load(savefile1)  ## should create Allrho_WBMatrix_list
      if (!exists("Allrho_WBMatrix_list")) {
        stop("Allrho_WBMatrix_list not found in ", savefile1)
      }

      if (length(Allrho_WBMatrix_list) < rr) {
        cat("       - Allrho_WBMatrix_list has only",
            length(Allrho_WBMatrix_list),
            "elements, but rr =", rr, "\n")
        rm(Allrho_WBMatrix_list)
        next
      }

      beta_chr <- Allrho_WBMatrix_list[[rr]]

      ## Your original structure: cbind t() across chromosomes
      beta_chr_t <- t(beta_chr)

      if (is.null(AllBeta)) {
        AllBeta <- beta_chr_t
      } else {
        AllBeta <- cbind(AllBeta, beta_chr_t)
      }

      rm(Allrho_WBMatrix_list, beta_chr, beta_chr_t)
      gc()
    }  ## end for chr

    if (is.null(AllBeta)) {
      cat("    - No Beta data for rho =", rho, " tauuse =", tauuse, "\n")
      Allrho_ZMatrix_list[[rr]] <- NULL
    } else {
      Allrho_ZMatrix_list[[rr]] <- AllBeta
    }

    rm(AllBeta)
    gc()
  }  ## end for rr (rho loop)

  ## (Optional) Save the combined Allrho_ZMatrix_list for this tauuse
  savefile_tau <- file.path(
    output_sub_folder_Beta,
    paste0("sim", sim,
           "_tau", tauuse,
           "_allchrs.RData")
  )
  save(Allrho_ZMatrix_list, file = savefile_tau)

  ## ----------------------------------
  ## 6b. For this tauuse, compute rowSums(log(.)) and track max
  ## ----------------------------------
  results_list <- vector("list", length(rho_vec))

  for (rr in seq_along(rho_vec)) {
    rho <- rho_vec[rr]
    combined_matrix <- Allrho_ZMatrix_list[[rr]]

    if (is.null(combined_matrix)) {
      cat("    - Skipping rho =", rho,
          " for tauuse =", tauuse, " (no combined matrix)\n")
      next
    }

    ## Avoid log(0) by adding epsilon where entries are 0
    adjusted_matrix <- ifelse(combined_matrix == 0,
                              combined_matrix + epsilon,
                              combined_matrix)

    sum_of_logs <- rowSums(log(adjusted_matrix))

    max_index  <- which.max(sum_of_logs)
    max_value  <- sum_of_logs[max_index]

    results_list[[rr]] <- list(
      Row_index           = max_index,
      Largest_sum_of_logs = max_value,
      Corresponding_rho   = rho,
      Corresponding_tauuse = tauuse
    )
  }

  ## Filter out any NULL entries (in case of missing data)
  non_empty <- !vapply(results_list, is.null, logical(1))
  if (!any(non_empty)) {
    cat(">>> No valid results for tauuse =", tauuse, "\n\n")
    rm(Allrho_ZMatrix_list, results_list)
    gc()
    next
  }

  results_list_ne <- results_list[non_empty]

  largest_values_tt  <- sapply(results_list_ne, `[[`, "Largest_sum_of_logs")
  largest_indices_tt <- sapply(results_list_ne, `[[`, "Row_index")
  largest_rhos_tt    <- sapply(results_list_ne, `[[`, "Corresponding_rho")
  largest_tau_tt     <- sapply(results_list_ne, `[[`, "Corresponding_tauuse")

  max_largest_value_tt <- max(largest_values_tt)
  max_index_tt         <- which(largest_values_tt == max_largest_value_tt)[1]

  corresponding_rho_tt    <- largest_rhos_tt[max_index_tt]
  corresponding_tauuse_tt <- largest_tau_tt[max_index_tt]
  corresponding_row_tt    <- largest_indices_tt[max_index_tt]

  ## Append to global trackers
  largest_values  <- c(largest_values,  max_largest_value_tt)
  largest_indices <- c(largest_indices, corresponding_row_tt)
  largest_rhos    <- c(largest_rhos,    corresponding_rho_tt)
  largest_tauuse  <- c(largest_tauuse,  corresponding_tauuse_tt)

  cat(">>> For tauuse =", tauuse,
      " best (within this tau) log-sum =", max_largest_value_tt,
      " at row =", corresponding_row_tt,
      " rho =", corresponding_rho_tt, "\n\n")

  rm(Allrho_ZMatrix_list, results_list, results_list_ne)
  gc()

}  ## end for tt over Tau columns

## ============================================================
## 7. Overall maximum across all tauuse & rho
## ============================================================
if (length(largest_values) == 0) {
  stop("No valid largest_values found")
}

overall_max_index      <- which.max(largest_values)
overall_largest_value  <- largest_values[overall_max_index]
overall_largest_index  <- largest_indices[overall_max_index]
overall_corresponding_rho <- largest_rhos[overall_max_index]
overall_tauuse         <- largest_tauuse[overall_max_index]

cat("=====================================================\n")
cat(">>> OVERALL MAX (3 traits, 1 pop: EUR)\n")
cat(">>> sim =", sim,
    " overall Largest_sum_of_logs =", overall_largest_value,
    " at row =", overall_largest_index,
    " rho =", overall_corresponding_rho,
    " tauuse =", overall_tauuse, "\n")
cat("=====================================================\n")

## ============================================================
## 8. Output as data.frame and write to CSV
## ============================================================
out_df <- data.frame(
  sim                 = sim,
  Largest_Sum_Of_Logs = overall_largest_value,
  Largest_Index       = overall_largest_index,
  Corresponding_Rho   = overall_corresponding_rho,
  Corresponding_Tauuse = overall_tauuse
)

## Save summary for this sim + scenario
summary_file <- file.path(
  sc_dir,
  paste0("MaxLogSum_3traits_sim", sim, "_scenario", scenarioIndex, ".csv")
)

write.table(out_df, file = summary_file,
            sep = ",", row.names = FALSE, quote = FALSE)

print(out_df)
