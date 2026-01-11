#!/usr/bin/env Rscript
rm(list=ls())

suppressPackageStartupMessages({
  library(MASS)
  library(mvtnorm)
  library(gtools)
})

# -----------------------------------------------------------
# Base directory (scenario folders live here)
# -----------------------------------------------------------
popvec <- c("AFR","EAS","EUR")

dirBase <- "/home/yatah3/projects/def-thchlava/yatah3/simulation/SimuGenotype/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0)
  stop("No scenario folders found under: ", dirBase)

# -----------------------------------------------------------
# Scenario grid (must match the simulation step)
# -----------------------------------------------------------
rho_grid <- c(0.2, 0.4, 0.6, 0.8, 0.9)
hsq_grid <- c(0.3, 0.5, 0.6, 0.8, 0.9)

train_grid <- list(
  c(5000,  5000,  50000),
  c(10000, 10000,  80000),
  c(15000, 15000, 100000),
  c(20000, 20000, 150000),
  c(25000, 25000, 200000)
)

val_grid <- list(
  c(5000,  5000,  5000),
  c(10000, 10000, 10000),
  c(15000, 15000, 15000),
  c(20000, 20000, 20000),
  c(25000, 25000, 25000)
)

# -----------------------------------------------------------
# helper to find scenario index from folder name (based on rho & hsq)
# -----------------------------------------------------------
get_scenario_index <- function(folder_name) {
  pattern <- "sim_hsq([0-9.]+)_rho([0-9.]+)"
  matches <- regmatches(folder_name, regexec(pattern, folder_name))[[1]]
  if (length(matches) == 3) {
    hsq_val <- as.numeric(matches[2])
    rho_val <- as.numeric(matches[3])
    idx <- which(abs(hsq_grid - hsq_val) < 1e-8 & abs(rho_grid - rho_val) < 1e-8)
    if (length(idx) == 1) return(idx)
  }
  return(NA_integer_)
}

# -----------------------------------------------------------
# main constants (same as your old script)
# -----------------------------------------------------------
nsim <- 20

# -----------------------------------------------------------
# MAIN loop over scenario folders
# -----------------------------------------------------------
for (sc_dir in scenario_dirs) {

  sc_name <- basename(sc_dir)
  idx <- get_scenario_index(sc_name)

  if (is.na(idx)) {
    warning("Could not match scenario index for folder: ", sc_name)
    next
  }

  valNsam_vec <- val_grid[[idx]]
  names(valNsam_vec) <- popvec

  cat("\n==============================\n")
  cat("Processing Scenario", idx, ":", sc_name, "\n")
  cat("Validation sizes (per pop): ",
      paste0(popvec, "=", valNsam_vec, collapse = ", "), "\n")
  cat("==============================\n")

  # ---------------------------------------------------------
  # Load true SNP info (contains hsq) - scenario-specific
  # ---------------------------------------------------------
  trueSNPfile <- paste0(sc_dir, "/trueSNP_setting.RData")
  if (!file.exists(trueSNPfile)) stop("Missing: ", trueSNPfile)

  load(file = trueSNPfile)  # must load hsq

  # ---------------------------------------------------------
  # Loop over populations EXACTLY like old script
  # ---------------------------------------------------------
  for (iiIndex in c(1:3)) {

    pop <- popvec[iiIndex]

    # In your old script, fam comes from dataname ".fam"
    dataname <- paste0(dirBase, pop, "AllChrs_bedformat")
    famfile  <- paste0(dataname, ".fam")
    if (!file.exists(famfile)) stop("Missing fam file: ", famfile)

    # In your scenario folder, we expect temp_<pop>.profile exists already
    tempfile <- paste0(sc_dir, "/temp_", pop)
    proffile <- paste0(tempfile, ".profile")
    if (!file.exists(proffile)) {
      warning("Skipping ", pop, " because profile missing: ", proffile)
      next
    }

    # OLD STYLE: loop sim=1..nsim
    for (sim in 1:nsim) {

      fam <- read.table(file = famfile, as.is = TRUE)

      xb <- y <- read.table(file = proffile, as.is = TRUE, skip = 1)[, c(1,2,6)]
      y[,3] <- y[,3]/sqrt(var(y[,3])) * sqrt(sum(hsq)) + rnorm(nrow(y))*sqrt(1-sum(hsq))

      print(cor(xb[,3], y[,3])^2)

      # ---- OLD FILENAMES (NO UNDERSCORES) ----
      testfile <- paste0(sc_dir, "/TestDid", pop, "sim", sim, ".txt")
      valfile  <- paste0(sc_dir, "/ValiDid", pop, "sim", sim, ".txt")
      disfile  <- paste0(sc_dir, "/Discoveryid", pop, "sim", sim, ".txt")
      tvfile   <- paste0(sc_dir, "/Test_ValiDid", pop, "sim", sim, ".txt")

      # ---- OLD deterministic split structure ----
      samids <- 1:nrow(fam)

      tidp <- samids[1:(valNsam_vec[iiIndex]/2)]
      vidp <- samids[((valNsam_vec[iiIndex]/2) + 1):valNsam_vec[iiIndex]]
      didp <- samids[(valNsam_vec[iiIndex] + 1):nrow(fam)]

      tlid  <- fam[tidp,1:2]
      valid <- fam[vidp,1:2]
      Disid <- fam[didp,1:2]

      testvalid <- fam[c(tidp, vidp), 1:2]

      write.table(tlid,      file=testfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
      write.table(valid,     file=valfile,  quote=FALSE, row.names=FALSE, col.names=FALSE)
      write.table(Disid,     file=disfile,  quote=FALSE, row.names=FALSE, col.names=FALSE)
      write.table(testvalid, file=tvfile,   quote=FALSE, row.names=FALSE, col.names=FALSE)

      # ---- OLD phenotype file names (NO UNDERSCORES) ----
      ytestvalfile <- paste0(sc_dir, "/test_validation.pheno_", pop, "sim", sim)
      ydisfile     <- paste0(sc_dir, "/training.pheno_", pop, "sim", sim)

      write.table(y[c(tidp,vidp),], file=ytestvalfile, row.names=FALSE, col.names=FALSE, quote=FALSE)
      write.table(y[didp,],         file=ydisfile,     row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
  }
}

cat("\nDone: OLD-style deterministic splits + OLD filenames, but scenario-aware folders.\n")
