rm(list = ls())

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
})

# =======================================================================
# Inputs: sim, scenarioIndex read from job array file
# =======================================================================
args  <- commandArgs(trailingOnly = TRUE)
jset  <- as.numeric(args[1])

file.rjobs <- "/lustre10/scratch/yatah3/yatah3/simulation/project1/step3/"
inputs     <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE)

scenarioIndex <- inputs[jset, 1]
sim           <- inputs[jset, 2]

cat(">>> Running ES-ALPHA extraction, sim =", sim, " scenario =", scenarioIndex, "\n")

# =======================================================================
# Directories
# =======================================================================
dirBase <- "/lustre10/scratch/yatah3/yatah3/simulation/SimuGenotype/three_traits/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)

cat(">>> Scenario:", sc_name, "\n")

# Beta matrices from step1
beta_dir <- file.path(sc_dir, "Trans2-new-3traits")

# Tau file
Taufile <- file.path(sc_dir, "EUR3traits_Tau_info_v2.RData")
load(Taufile)   # loads AbsTauvec

# =======================================================================
# Fixed parameters used in Step 1 (must match Step 1 output filenames)
# =======================================================================
penalty     <- "RealmixLOG"
warmStart   <- 1
singleStart <- 1
Zscale      <- 1

# =======================================================================
# Load required functions
# =======================================================================
functions_folder <- "/lustre10/scratch/yatah3/yatah3/simulation/Rfunction/"
source(paste0(functions_folder, "AllPRS_Rfunctions.r"))
source(paste0(functions_folder, "PRS_utility.r"))
source(paste0(functions_folder, "Iterative_Rfunctions.r"))
source(paste0(functions_folder, "PlinkLD_transform.R"))

# =======================================================================
# Read MaxLogSum_* to get the correct tauuse + row index
# =======================================================================
summary_file <- file.path(
  sc_dir,
  paste0("MaxLogSum_3traits_sim", sim, "_scenario", scenarioIndex, ".csv")
)

best_res <- read.csv(summary_file)

tauuse       <- best_res$Corresponding_Tauuse
row_index    <- best_res$Largest_Index

cat(">>> Using optimal tauuse =", tauuse, 
    "optimal alpha index =", row_index, "\n")

# =======================================================================
# Load Beta matrices for ALL chromosomes using this tauuse
# =======================================================================

for (chr in 1:22) {

  savefile1 <- paste0(
    beta_dir, "/", penalty,
    "chr", chr, "_3traits",
    "warmStart", warmStart,
    "sim", sim,
    "Zscale", Zscale,
    "tauuse", tauuse,
    "singleStart", singleStart,
    ".RData"
  )

  if (!file.exists(savefile1)) {
    cat(">>> No file for chr", chr, "skipping \n")
    next
  }

  BetaMatrixS <- Gen_One_BetaMatrix(savefile1, 3, 1)[[1]]
  BetaMatrix3 <- Gen_One_BetaMatrix(savefile1, 3, 2)[[1]]
  BetaMatrixE <- Gen_One_BetaMatrix(savefile1, 3, 3)[[1]]

 
  if (chr == 1) {
    AllBetaMatrixS <- BetaMatrixS
    AllBetaMatrix3 <- BetaMatrix3
    AllBetaMatrixE <- BetaMatrixE
  } else {
    AllBetaMatrixS <- cbind(AllBetaMatrixS, BetaMatrixS)
    AllBetaMatrix3 <- cbind(AllBetaMatrix3, BetaMatrix3)
    AllBetaMatrixE <- cbind(AllBetaMatrixE, BetaMatrixE)
  }

  rm(BetaMatrixS, BetaMatrix3, BetaMatrixE)
  gc()
}

BetaS_all <- AllBetaMatrixS
Beta3_all <- AllBetaMatrix3
BetaE_all <- AllBetaMatrixE


total_snps <- nrow(BetaS_all)
cat(">>> Total SNPs =", total_snps, "\n")

# =======================================================================
# Compute ES-ALPHA matrix
# =======================================================================
es_alphaMatrix <- matrix(0, total_snps, 8)

for (i in 1:total_snps) {

  v1 <- BetaS_all[i,]
  v2 <- Beta3_all[i,]
  v3 <- BetaE_all[i,]

  pi_0   <- mean(v1==0 & v2==0 & v3==0)
  pi_1   <- mean(v1!=0 & v2==0 & v3==0)
  pi_2   <- mean(v1==0 & v2!=0 & v3==0)
  pi_3   <- mean(v1==0 & v2==0 & v3!=0)
  pi_12  <- mean(v1!=0 & v2!=0 & v3==0)
  pi_13  <- mean(v1!=0 & v2==0 & v3!=0)
  pi_23  <- mean(v1==0 & v2!=0 & v3!=0)
  pi_123 <- mean(v1!=0 & v2!=0 & v3!=0)

  es_alphaMatrix[i,] <- c(pi_0,pi_1,pi_2,pi_3,pi_12,pi_13,pi_23,pi_123)
}

# =======================================================================
# Extract selected row
# =======================================================================
selected_row <- es_alphaMatrix[row_index, ]
cat(">>> Selected alpha row:\n")
print(selected_row)

# =======================================================================
# Save results
# =======================================================================
outfile <- file.path(
  sc_dir,
  paste0("es_alpha_best_sim", sim,
         "_scenario", scenarioIndex,
         "_row", row_index,
         "_tau", tauuse, ".RData")
)

save(es_alphaMatrix, row_index, tauuse, selected_row, file = outfile)

cat(">>> Saved final output:", outfile, "\n")
