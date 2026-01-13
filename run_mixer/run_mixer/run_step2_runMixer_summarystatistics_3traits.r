#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(SummaryLasso)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(jsonlite)  # for fromJSON()
})

# ============================================================
# 1) Job inputs (from SLURM array)
#    Same format as previous MiXeR step:
#    input1.txt:  scenarioIndex  sim
# ============================================================
args  <- commandArgs(trailingOnly = TRUE)
jset <- as.numeric(args[1])

file.rjobs <- "/lustre06/project/6005709/yatah3/simulation/run_mixer/"
inputs <- read.table(paste0(file.rjobs, "input1.txt"),
                     header = TRUE, as.is = TRUE)

scenarioIndex <- inputs[jset, 1]
sim           <- inputs[jset, 2]

cat(">>> sim =", sim, " scenario =", scenarioIndex, "\n")

# one population, 3 traits (same as previous step)
popvec       <- "EUR"
trait_names  <- c("trait1", "trait2", "trait3")

# ============================================================
# 2) Scenario selection (same dir structure as previous step)
# ============================================================
dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/three_traits/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0)
  stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs))
  stop("scenarioIndex exceeds available scenarios (", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)
cat(">>> Using scenario:", sc_name, "\n")

# ------------------------------------------------------------
# 3) Paths to MiXeR output (input) and summary CSV (output)
# ------------------------------------------------------------
# MiXeR step created this folder:
#   sc_dir / Mixer_sim{sim}
output_sub_folder <- file.path(sc_dir, paste0("Mixer_sim", sim))

# Post-MiXeR summaries will go to:
#   sc_dir / outputvisualonepop_sim{sim} /
output_sub_folder2 <- file.path(sc_dir, paste0("outputvisualonepop_sim", sim))
dir.create(output_sub_folder2, recursive = TRUE, showWarnings = FALSE)

cat("MiXeR output folder:        ", output_sub_folder,  "\n")
cat("Post-MiXeR summary folder:  ", output_sub_folder2, "\n")

# ------------------------------------------------------------
# 4) Read MiXeR .fit.rep.json for each trait and write CSV
# ------------------------------------------------------------
for (ii in 1:3) {

  # Example filename: trait1_sim9.fit.rep.json
  inputfile3 <- file.path(
    output_sub_folder,
    paste0(trait_names[ii], "_sim", sim, ".fit.rep.json")
  )

  # Example output: trait1_sim9.csv
  outfile3 <- file.path(
    output_sub_folder2,
    paste0(trait_names[ii], "_sim", sim, ".csv")
  )

  if (!file.exists(inputfile3)) {
    cat("  [WARN] Missing file:", inputfile3, "\n")
    next
  }

  cat("  Reading:", inputfile3, "\n")
  jdata <- fromJSON(inputfile3)

  # Safely extract parameters (set NA if missing)
  pi_val        <- if (!is.null(jdata$params$pi))        jdata$params$pi        else NA
  sig2_beta_val <- if (!is.null(jdata$params$sig2_beta)) jdata$params$sig2_beta else NA
  sig2_zero_val <- if (!is.null(jdata$params$sig2_zero)) jdata$params$sig2_zero else NA

  # Clean one-row data frame
  result <- data.frame(
    Population = popvec,
    Simulation = sim,
    Trait      = trait_names[ii],
    pi         = pi_val,
    sig2_beta  = sig2_beta_val,
    sig2_zero  = sig2_zero_val,
    stringsAsFactors = FALSE
  )

  write.csv(result, outfile3, row.names = FALSE)
  cat("  [OK] Saved:", outfile3, "\n")
}

cat("\nPost-MiXeR summary finished for sim", sim,
    "scenario", scenarioIndex, "\n")
