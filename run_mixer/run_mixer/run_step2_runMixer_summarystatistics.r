rm(list = ls())
library(dplyr)
library(stringr)

library(jsonlite)


# -------------------------------------------------------------------
# 1. Parse job index (array ID)
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("❌ Missing argument: SLURM_ARRAY_TASK_ID")
jset <- as.numeric(args[1])

# -------------------------------------------------------------------
# 2. Load input table (scenarioIndex + sim)
# -------------------------------------------------------------------
file.rjobs <- "/lustre06/project/6005709/yatah3/simulation/run_mixer/"
inputs <- read.table(paste0(file.rjobs, "input1.txt"), header = TRUE, as.is = TRUE)


scenarioIndex <- inputs[jset, 1]
sim           <- inputs[jset, 2]

cat(" scenarioIndex =", scenarioIndex,
    " sim =", sim, "\n")

# -------------------------------------------------------------------
# 3. Locate correct scenario folder
# -------------------------------------------------------------------
dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("scenario", basename(scenario_dirs)) |
                               grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0)
  stop("No scenario folders found under:", dirBase)
if (scenarioIndex < 1 || scenarioIndex > length(scenario_dirs))
  stop("Invalid scenarioIndex:", scenarioIndex, "(must be between 1 and", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)
cat(" Scenario folder selected:", sc_name, "\n")

# -------------------------------------------------------------------
# 4. Define constants and paths
# -------------------------------------------------------------------
dirOutput <- sc_dir
dir1000G  <- paste0(dirOutput, "/ld_mixer/")
output_sub_folder <- file.path(sc_dir, paste0("Mixer_sim", sim))
system(paste0("mkdir -p ", output_sub_folder))

# Mixer python + lib paths
dirPython <- "/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/precimed/mixer.py"
lib       <- "/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/src/build/lib/libbgmg.so"
popvec <- c("AFR", "EAS", "EUR")

# -------------------------------------------------------------------
# Post-MiXeR: compute summary statistics for this sim
# -------------------------------------------------------------------
output_sub_folder2 <- paste0(sc_dir, "/outputvisualonepop_sim", sim, "/")
system(paste0("mkdir -p ", output_sub_folder2))

for (ii in 1:3) {
  inputfile3 <- paste0(output_sub_folder, "/", popvec[ii], "_sim", sim, ".fit.rep.json")
  outfile3   <- paste0(output_sub_folder2, "/", popvec[ii], "_sim", sim, ".csv")

  if (!file.exists(inputfile3)) {
    cat("  ⚠️ Missing file:", inputfile3, "\n")
    next
  }

  cat("  Reading:", inputfile3, "\n")
  jdata <- fromJSON(inputfile3)

  # Safely extract parameters (set NA if missing)
  pi_val        <- if (!is.null(jdata$params$pi)) jdata$params$pi else NA
  sig2_beta_val <- if (!is.null(jdata$params$sig2_beta)) jdata$params$sig2_beta else NA
  sig2_zero_val <- if (!is.null(jdata$params$sig2_zero)) jdata$params$sig2_zero else NA

  # Build a clean one-row data frame
  result <- data.frame(
    Population = popvec[ii],
    Simulation = sim,
    pi = pi_val,
    sig2_beta = sig2_beta_val,
    sig2_zero = sig2_zero_val
  )

  # Save as CSV (comma-separated)
  write.csv(result, outfile3, row.names = FALSE)
  cat("  ✅ Saved:", outfile3, "\n")
}
