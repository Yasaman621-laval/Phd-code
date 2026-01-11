#!/usr/bin/env Rscript
rm(list = ls())

# ============================================================
# 0. Libraries
# ============================================================
library(SummaryLasso)
library(MASS)
library(mvtnorm)
library(gtools)
library(dplyr)
library(data.table)
library(stringr)

# ============================================================
# 1. Parse job index (SLURM array ID)
#    sim, runIndex, ttIndex, rho, scenarioIndex
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing argument: SLURM_ARRAY_TASK_ID")

jset <- as.numeric(args[1])

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project2/step3/"
inputs     <- read.table(paste0(file.rjobs, "input.txt"),
                         as.is = TRUE, header = TRUE)

sim           <- inputs[jset, 1]
runIndex      <- inputs[jset, 2]
ttIndex       <- inputs[jset, 3]
rr           <- inputs[jset, 4]
scenarioIndex <- inputs[jset, 5]

cat(">>> sim =", sim,
    " runIndex =", runIndex,
    " ttIndex =", ttIndex,
    " rr =", rr,
    " scenario =", scenarioIndex, "\n")

# ============================================================
# 2. Scenario detection
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0)
  stop("No scenario directories found under: ", dirBase)

if (scenarioIndex > length(scenario_dirs))
  stop("scenarioIndex exceeds available scenarios")

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)

cat(">>> Using scenario:", sc_name, "\n")

# ============================================================
# 3. Training sample sizes

# ============================================================
train_str       <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)

TrainingNsam <- as.numeric(
  sapply(strsplit(train_str_clean, "-")[[1]],
         function(x) eval(parse(text = x)))
)

cat(">>> TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")

# ============================================================
# 4. General settings
# ============================================================
Zscale   <- 1
nsim     <- 10
savename <- "R2output"
penalty  <- "RealmixLOG"

dirtemp      <- Sys.getenv("SLURM_TMPDIR")
folder_predi <- paste0(dirtemp, "/")

ordersequse_vec <- c(1, 1)
all_usedtrait   <- c("1,3", "2,3")


functionsfolder <- "/lustre09/project/6005709/yatah3/simulation/Rfunction/"


dirSimuData     <- paste0(sc_dir, "/")


dirSimuInOutput <- sc_dir

output_sub_folder <- file.path(dirSimuInOutput, savename, "/")
system(paste("mkdir -p", output_sub_folder))


output_sub_folder_uDensity <- file.path(dirSimuInOutput, "new_densityU/")


R2output <- paste0(output_sub_folder, "R2output1/")
system(paste("mkdir -p", R2output))


Scoreoutput <- paste0(folder_predi, "/Scoreoutput1/")
system(paste("mkdir -p", Scoreoutput))

# ============================================================
# 5. Load helper functions
# ============================================================
source(paste0(functionsfolder, "AllPRS_Rfunctions.r"))
source(paste0(functionsfolder, "PRS_utility.r"))
source(paste0(functionsfolder, "Iterative_Rfunctions.r"))
source(paste0(functionsfolder, "PlinkLD_transform.R"))

# ============================================================
# 6. Basic setup
# ============================================================
popvec <- c("AFR", "EAS", "EUR")

warmStart <- 1

refAllelefile <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/refereceAlleleSNPs.txt"
refAllelefile_list <- list()
refAllelefile_list[[1]] <- refAllelefile
refAllelefile_list[[2]] <- refAllelefile

K <- 2

ordersequse   <- ordersequse_vec[runIndex]
usedtrait     <- all_usedtrait[runIndex]
usedtraitsvec <- unlist(strsplit(as.character(usedtrait), split = ","))
usedtraitIndex <- as.numeric(usedtraitsvec)

mainindex <- iiIndex0 <- usedtraitIndex[1]
iiIndexT  <- usedtraitIndex[2]

phenoIndex <- 3
gindex     <- 1
plinkver   <- 2
cindex     <- 1
Before     <- 0

popuseY <- popvec[mainindex]

cindex  <- 1
N       <- TrainingNsam[3]   # same logic as original
numsnp <- 400

rho_vec <- c(seq(0, 0.9, 0.1), 0.95)

singleStart <- 1
PV          <- 2

Taufile <- paste0(sc_dir, "/", popvec[mainindex], "_",
                  popvec[iiIndexT], "Tau_info_v2.RData")
if (!file.exists(Taufile))
  stop("Taufile missing: ", Taufile)

load(Taufile)   # loads AbsTauvec




# ============================================================
# 7. Main loop (ttIndex only, rho fixed via jset)
# ============================================================
rho <- rho_vec[rr]

tt = ttIndex

  tauuse <- AbsTauvec[sim, tt]
  cat(">>> tau column =", tt, " -> tauuse =", tauuse, "\n")



  cat(">>> rho =", rho, "\n")

  savecore <- paste0("usedtrait", usedtrait,
                     "sim", sim,
                     "tauuse", tauuse,
                     "rho", rho)

  savef2_1 <- paste0(R2output, penalty,

                     "usedtrait", usedtrait,
                     "warmStart", warmStart,
                     "sim", sim,
                     "Zscale", Zscale,
                     "singleStart", singleStart,
                     "tauuse", tauuse,
                     "rho", rho,
                     "_wB1R2.RData")

  savef2_2 <- paste0(R2output, penalty,

                     "usedtrait", usedtrait,
                     "warmStart", warmStart,
                     "sim", sim,
                     "Zscale", Zscale,
                     "singleStart", singleStart,
                     "tauuse", tauuse,
                     "rho", rho,
                     "_wB2R2.RData")

    savef2_list <- list()
    savef2_list[[1]] <- savef2_1
    savef2_list[[2]] <- savef2_2

  output_sub_folder_predi <- paste0(folder_predi, savecore, "/")
  system(paste("mkdir -p", output_sub_folder_predi))

  # ----------------------------------------------------------
  # Collect AllBetaMatrix across chromosomes
  # ----------------------------------------------------------
  for (chr in 1:22) {

      subtau_saveoutfile <- paste0(output_sub_folder_uDensity, penalty, "_chr", chr, 
                                 "_usedtrait", usedtrait, "_warmStart", warmStart, 
                                 "_sim", sim, "_Zscale", Zscale, "_singleStart", singleStart, 
                                 "_tauuse", tauuse, "_DensityU.RData")


    if (!file.exists(subtau_saveoutfile)) {
      stop("Missing file: ", subtau_saveoutfile)

    }

    load(subtau_saveoutfile)  # loads Allrho_WBMatrix_list

    if (chr == 1) {
      AllBetaMatrix1 <- t(Allrho_WBMatrix_list[[rr]][[1]])
      AllBetaMatrix2 <- t(Allrho_WBMatrix_list[[rr]][[2]])
    } else {
      AllBetaMatrix1 <- cbind(AllBetaMatrix1,
                              t(Allrho_WBMatrix_list[[rr]][[1]]))
      AllBetaMatrix2 <- cbind(AllBetaMatrix2,
                              t(Allrho_WBMatrix_list[[rr]][[2]]))
    }

    rm(Allrho_WBMatrix_list)
    gc()
  }

  if (!identical(colnames(AllBetaMatrix1), colnames(AllBetaMatrix2)))
    stop("Column mismatch between Beta matrices")



  if (anyNA(AllBetaMatrix1) || anyNA(AllBetaMatrix2))
    stop("NA detected in Beta matrices")

    AllBetaMatrix_list <- list()
    AllBetaMatrix_list[[1]] <- AllBetaMatrix1
    AllBetaMatrix_list[[2]] <- AllBetaMatrix2
  rm(AllBetaMatrix1, AllBetaMatrix2)
  gc()


  GenPreR2_Chrs_saveScore(
    AllBetaMatrix_list,
    savef2_list,
    refAllelefile_list,
    output_sub_folder_predi,
    PV,
    popuseY,
    dirSimuData,
    sim,
    savename,
    Scoreoutput,
    savecore
  )

  rm(AllBetaMatrix_list)
  gc()


cat("Completed sim", sim,
    "runIndex", runIndex,
    "ttIndex", ttIndex,
    "rho", rho,
    "scenario", scenarioIndex, "\n")
