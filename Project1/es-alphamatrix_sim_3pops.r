#!/usr/bin/env Rscript

rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# =======================================================================
# 0) Inputs: sim, scenarioIndex read from job array file
# =======================================================================
args  <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("? Missing argument: SLURM_ARRAY_TASK_ID")
jset  <- as.numeric(args[1])

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project1/step3/"
inputs     <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE)

# Expected columns: scenarioIndex, sim
scenarioIndex <- inputs[jset, 1]
sim           <- inputs[jset, 2]

cat(">>> Running ES-ALPHA extraction (3-pop) sim =", sim,
    " scenarioIndex =", scenarioIndex, "\n")

# =======================================================================
# 1) Scenario detection (must match 3-pop MaxLogSum step)
# =======================================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"
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

cat(">>> Scenario directory:", sc_dir, "\n")

# =======================================================================
# 2) Fixed parameters (must match the step that created Beta matrices)
# =======================================================================
penalty     <- "RealmixLOG"
warmStart   <- 1
singleStart <- 1
Zscale      <- 1

# 3-population order (same as in MaxLogSum step)
popvec <- c("AFR", "EAS", "EUR")
K      <- length(popvec)

# Directory where Beta RData lives (EDIT if your folder name differs)
# This should match the folder used when you ran transConditionalU
# to generate the Beta matrices.
beta_dir <- file.path(sc_dir, "Trans2-new-3")

cat(">>> Beta directory:", beta_dir, "\n")

# =======================================================================
# 3) Load functions (if Gen_One_BetaMatrix is defined there)
# =======================================================================
functions_folder <- "/lustre09/project/6005709/yatah3/simulation/Rfunction/"
source(paste0(functions_folder, "AllPRS_Rfunctions.r"))
source(paste0(functions_folder, "PRS_utility.r"))
source(paste0(functions_folder, "Iterative_Rfunctions.r"))
source(paste0(functions_folder, "PlinkLD_transform.R"))

# =======================================================================
# 4) Read MaxLogSum_sim* to get optimal tauuse + alpha row index
# =======================================================================
summary_file <- file.path(
  sc_dir,
  paste0("MaxLogSum_sim", sim, "_scenario", scenarioIndex, ".csv")
)

if (!file.exists(summary_file)) {
  stop("Cannot find MaxLogSum summary file: ", summary_file)
}

best_res <- read.csv(summary_file, header = TRUE)

tauuse    <- best_res$Corresponding_Tauuse[1]
row_index <- best_res$Largest_Index[1]

cat(">>> Using optimal tauuse =", tauuse,
    "optimal alpha index =", row_index, "\n")

# =======================================================================
# 5) Load Beta matrices for ALL chromosomes at this tauuse
# =======================================================================

for (chr in 1:22) {

  # This naming must match the file created in the 3-pop step1 run
  savefile1 <- file.path(
    beta_dir,
    paste0(
     penalty,"chr",chr,"_3pop","warmStart",warmStart,"sim",sim,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
  )

  if (!file.exists(savefile1)) {
    cat("No Beta file for chr", chr, " (", savefile1, ") skipping \n")
    next
  }

  # Gen_One_BetaMatrix(savefile1, K, kpop) should return list, [[1]] is matrix
  BetaAFR <- Gen_One_BetaMatrix(savefile1, K, 1)[[1]]
  BetaEAS <- Gen_One_BetaMatrix(savefile1, K, 2)[[1]]
  BetaEUR <- Gen_One_BetaMatrix(savefile1, K, 3)[[1]]

 
  if (chr == 1) {
    AllBetaAFR <- BetaAFR
    AllBetaEAS <- BetaEAS
    AllBetaEUR <- BetaEUR
  } else {
    AllBetaAFR <- cbind(AllBetaAFR, BetaAFR)
    AllBetaEAS <- cbind(AllBetaEAS, BetaEAS)
    AllBetaEUR <- cbind(AllBetaEUR, BetaEUR)
  }

  rm(BetaAFR, BetaEAS, BetaEUR)
  gc()
}

BetaAFR_all <- AllBetaAFR
BetaEAS_all <- AllBetaEAS
BetaEUR_all <- AllBetaEUR


# number of "alpha-samples" (rows) is the same for all three matrices
n_alpha <- nrow(BetaAFR_all)

cat(">>> Combined Beta dimensions:\n")
cat("    AFR:", dim(BetaAFR_all), "\n")
cat("    EAS:", dim(BetaEAS_all), "\n")
cat("    EUR:", dim(BetaEUR_all), "\n")

# =======================================================================
# 6) Build ES-ALPHA matrix (pi_0,...,pi_123) across all SNPs
#     - each row i corresponds to one alpha index
#     - each column corresponds to pattern:
#       (0,0,0), (1,0,0), (0,1,0), (0,0,1),
#       (1,1,0), (1,0,1), (0,1,1), (1,1,1)
# =======================================================================
es_alphaMatrix <- matrix(NA_real_, nrow = n_alpha, ncol = 8)
colnames(es_alphaMatrix) <- c(
  "pi_0",  "pi_1",  "pi_2",  "pi_3",
  "pi_12", "pi_13", "pi_23", "pi_123"
)

for (i in seq_len(n_alpha)) {

  v1 <- BetaAFR_all[i, ]  # AFR
  v2 <- BetaEAS_all[i, ]  # EAS
  v3 <- BetaEUR_all[i, ]  # EUR

  pi_0   <- mean(v1 == 0 & v2 == 0 & v3 == 0)
  pi_1   <- mean(v1 != 0 & v2 == 0 & v3 == 0)
  pi_2   <- mean(v1 == 0 & v2 != 0 & v3 == 0)
  pi_3   <- mean(v1 == 0 & v2 == 0 & v3 != 0)
  pi_12  <- mean(v1 != 0 & v2 != 0 & v3 == 0)
  pi_13  <- mean(v1 != 0 & v2 == 0 & v3 != 0)
  pi_23  <- mean(v1 == 0 & v2 != 0 & v3 != 0)
  pi_123 <- mean(v1 != 0 & v2 != 0 & v3 != 0)

  es_alphaMatrix[i, ] <- c(
    pi_0, pi_1, pi_2, pi_3,
    pi_12, pi_13, pi_23, pi_123
  )
}

# =======================================================================
# 7) Extract the selected row (the optimal alpha index)
# =======================================================================
if (row_index < 1 || row_index > n_alpha) {
  stop("row_index out of bounds: ", row_index,
       " (n_alpha = ", n_alpha, ")")
}

selected_row <- es_alphaMatrix[row_index, ]

cat(">>> Selected alpha row (index", row_index, "):\n")
print(selected_row)

# =======================================================================
# 8) Save results
# =======================================================================
outfile <- file.path(
  sc_dir,
  paste0("es_alpha_best_sim", sim,
         "_scenario", scenarioIndex,
         "_row", row_index,
         "_tau", tauuse, ".RData")
)

save(es_alphaMatrix, row_index, tauuse, selected_row, file = outfile)

cat(">>> Saved final ES-ALPHA output:", outfile, "\n")
