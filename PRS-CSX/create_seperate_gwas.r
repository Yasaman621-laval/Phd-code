#!/usr/bin/env Rscript
rm(list = ls())
library(data.table)
library(dplyr)
library(tidyr)

cat("\n===== PRS-CSx SUMSTATS BUILDER STARTED =====\n")

# ---------------------------------------------------------
# 1. SLURM INPUT
# ---------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing SLURM_ARRAY_TASK_ID")
jset <- as.numeric(args[1])

input_file <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/input_sim_scenario.txt"
inputs <- read.table(input_file, header = TRUE, as.is = TRUE)

sim           <- inputs[jset, 2]
scenarioIndex <- inputs[jset, 1]

# ---------------------------------------------------------
# 2. Locate scenario directory
# ---------------------------------------------------------
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype"

scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (scenarioIndex > length(scenario_dirs))
  stop("scenarioIndex exceeds available simulation scenarios!")

sc_dir <- scenario_dirs[scenarioIndex]
cat(">>> Using simulation scenario:", sc_dir, "\n")

# ---------------------------------------------------------
# 3. Load multi-pop GWASbetaStandard file
# ---------------------------------------------------------
GWASbetafile <- file.path(sc_dir, paste0("GWASbetaStandard_allchrssim", sim, ".txt"))
cat("Reading:", GWASbetafile, "\n")

GWASbeta <- fread(GWASbetafile)
if (nrow(GWASbeta) == 0) stop("GWASbeta file is empty!")

# Expected columns (from your examples):
# CHR, SNP, A1, b1, SE1, p1, b2, SE2, p2, b3, SE3, p3, ...

# ---------------------------------------------------------
# 4. Derive A1 and A2 for PRS-CSx sumstats
# ---------------------------------------------------------
# SNP format: "1:752721_A_G"
parts   <- tstrsplit(GWASbeta$SNP, "_")
allele1 <- parts[[2]]
allele2 <- parts[[3]]

# A1 in GWASbeta is the effect allele (from your example)
A1_sum <- GWASbeta$A1

# A2 is "the other allele" in the SNP suffix
A2_sum <- ifelse(
  A1_sum == allele1, allele2,
  ifelse(A1_sum == allele2, allele1, NA_character_)
)

# Safety check: if any NA, something is inconsistent
if (any(is.na(A2_sum))) {
  warning("Some A2 could not be matched from SNP suffix; check GWASbeta A1 vs SNP pattern.")
}

# ---------------------------------------------------------
# 5. Helper to write one PRS-CSx-format sumstats file
# ---------------------------------------------------------
write_sumstats <- function(beta, p, outfile) {
  df <- data.frame(
    SNP  = GWASbeta$SNP,
    A1   = A1_sum,
    A2   = A2_sum,
    BETA = beta,
    P    = p
  )
  fwrite(df, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  wrote:", outfile, "\n")
}

# ---------------------------------------------------------
# 6. Write AFR, EAS, EUR sumstats
# ---------------------------------------------------------
write_sumstats(GWASbeta$b1, GWASbeta$p1,
               file.path(sc_dir, paste0("GWAS_AFR_", sim, ".tsv")))

write_sumstats(GWASbeta$b2, GWASbeta$p2,
               file.path(sc_dir, paste0("GWAS_EAS_", sim, ".tsv")))

write_sumstats(GWASbeta$b3, GWASbeta$p3,
               file.path(sc_dir, paste0("GWAS_EUR_", sim, ".tsv")))

cat("\n===== DONE: PRS-CSx SUMSTATS CREATED SUCCESSFULLY =====\n")
