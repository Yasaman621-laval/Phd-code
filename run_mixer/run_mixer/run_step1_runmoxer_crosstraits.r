#!/usr/bin/env Rscript
rm(list = ls())

suppressPackageStartupMessages({
  library(SummaryLasso)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ============================================================
# 1) Job inputs (from SLURM array)
# ============================================================
args  <- commandArgs(trailingOnly = TRUE)
jset <- as.numeric(args[1])

file.rjobs <- "/lustre06/project/6005709/yatah3/simulation/run_mixer/"
inputs <- read.table(paste0(file.rjobs, "input1.txt"),
                     header = TRUE, as.is = TRUE)

scenarioIndex <- inputs[jset, 1]
sim           <- inputs[jset, 2]

cat(">>> sim =", sim, " scenario =", scenarioIndex, "\n")

# ============================================================
# 2) Global config
# ============================================================
nsim       <- 10
savename   <- "Trans2-new-3traits"   # not really used here, but kept
penalty    <- "RealmixLOG"
NumL       <- 10
subNumL    <- 10
WeightN    <- 0
Zscale     <- 1
warmStart  <- 1
singleStart<- 1
NumIter    <- 1000

# one population here (EUR), 3 traits
popvec       <- "EUR"
trait_names  <- c("trait1", "trait2", "trait3")

# ============================================================
# 3) Scenario selection + TrainingNsam parsing
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

# output folder per scenario + sim
output_sub_folder <- file.path(sc_dir, paste0("Mixer_sim", sim))
dir.create(output_sub_folder, recursive = TRUE, showWarnings = FALSE)

# parse TrainingNsam from folder name suffix after "_train"
# example: "sim_hsq0.3_rho0.2_trainEUR50000"
train_str <- sub(".*_train", "", sc_name)      # "EUR50000"
train_num <- as.numeric(gsub("[^0-9]", "", train_str))
if (is.na(train_num)) stop("Could not parse TrainingNsam from scenario name: ", sc_name)

# *** IMPORTANT: 3 traits -> vector of length 3 ***
TrainingNsam <- rep(train_num, 3)

cat("    TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")

# -------------------------------------------------------------------
# 5. Create summary statistics for this simulation only
# -------------------------------------------------------------------
GWASbetafile <- file.path(sc_dir,
                          paste0("GWASbetaStandard_allchrs_3traits_sim", sim, ".txt"))
if (!file.exists(GWASbetafile))
  stop("Missing GWAS file: ", GWASbetafile)
cat(" Reading:", GWASbetafile, "\n")

GWASbeta <- read.delim(GWASbetafile, header = TRUE, sep = "\t", quote = "",
                       fill = TRUE, strip.white = TRUE, comment.char = "",
                       stringsAsFactors = FALSE)

# Add A2 if missing (parsed from SNP like "CHR:BP_A1_A2")
if (!"A2" %in% names(GWASbeta)) {
  GWASbeta <- GWASbeta %>%
    mutate(
      A2 = str_extract(SNP, "(?<=_[A-Za-z])_[A-Za-z]+$") %>%
           str_remove("^_")
    ) %>%
    select(
      CHR, SNP, A1, A2,
      Zobs1, b1, SE1, p1,
      Zobs2, b2, SE2, p2,
      Zobs3, b3, SE3, p3
    )
}

# Split per trait and write out per-sim files
for (ii in 1:3) {
  cols_needed <- c("CHR", "SNP", "A1",
                   paste0("b", ii),
                   paste0("SE", ii),
                   paste0("p", ii))
  missing_cols <- setdiff(cols_needed, names(GWASbeta))
  if (length(missing_cols) > 0) {
    cat(" Missing columns for trait", ii, ":", paste(missing_cols, collapse = ", "), "\n")
    next
  }

  tmp <- GWASbeta[, cols_needed]
  names(tmp) <- c("CHR", "SNP", "A1", "BETA", "SE", "P")
  outfile <- file.path(sc_dir,
                       paste0("summary_stats_trait", ii, "_sim", sim, ".txt"))
  write.table(tmp, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Created:", outfile, "\n")
}

# -------------------------------------------------------------------
# 6. Harmonize per-sim summary files
# -------------------------------------------------------------------
for (ii in 1:3) {
  infile  <- file.path(sc_dir,
                       paste0("summary_stats_trait", ii, "_sim", sim, ".txt"))
  outfile <- file.path(sc_dir,
                       paste0("summary_stats_trait", ii, "_sim", sim, "_harmonized.txt"))

  if (!file.exists(infile)) {
    cat("  [WARN] Missing input for harmonization:", infile, "\n")
    next
  }

  df <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
  names(df) <- toupper(names(df))

  df <- df %>%
    mutate(
      CHR = if ("CHR" %in% names(.)) CHR else as.numeric(sub(":.*", "", SNP)),
      BP  = as.numeric(sub(".*:(\\d+)_.*", "\\1", SNP)),
      A1  = toupper(A1),
      A2  = if ("A2" %in% names(.)) toupper(A2) else
              sub(".*_[A-Z]+_([A-Z]+)$", "\\1", SNP),
      Z   = ifelse(SE != 0, BETA / SE, 0),
      N   = TrainingNsam[ii]
    ) %>%
    select(SNP, CHR, BP, A1, A2, N, Z) %>%
    arrange(CHR, BP)

  if (any(is.na(df$N))) {
    cat("  [WARN] N is NA for trait", ii, "check TrainingNsam vector.\n")
  }

  write.table(df, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  cat("  Harmonized:", outfile, "\n")
}

# -------------------------------------------------------------------
# 7. Run MiXeR per trait for this sim
# -------------------------------------------------------------------
dir1000G  <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/ld_mixer/" # fixed
bim_dir   <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/sim_hsq0.5_rho0.4_train10000-10000-80000/"        # fixed

dirPython <- "/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/precimed/mixer.py"
lib       <- "/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/src/build/lib/libbgmg.so"

for (ii in 1:3) {
  cat("\n===========================================================\n")
  cat(">>> Processing trait:", trait_names[ii], "for sim", sim, "\n")
  cat("===========================================================\n")

  sumstatfile1 <- file.path(sc_dir,
                            paste0("summary_stats_trait", ii, "_sim", sim, "_harmonized.txt"))
  if (!file.exists(sumstatfile1)) {
    cat("  [WARN] Harmonized sumstats not found, skipping:", sumstatfile1, "\n")
    next
  }

  outfile <- file.path(output_sub_folder,
                       paste0(trait_names[ii], "_sim", sim, ".fit.rep"))

  # Paths for MiXeR inputs (fixed LD reference)
  snpfile <- paste0(dir1000G, popvec,
                    "_allchr.prune_maf0p05_rand2M_r2p8.snps")
  bimfile <- paste0(bim_dir,
                    "/Chr@/PlinkFormat/", popvec, "chr@bedformat.bim")
  ldfile  <- paste0(dir1000G,
                    "chr@.", popvec, ".ld")

  # ------------------------------------------------------------
  # OPTIONAL: quick overlap check to avoid arg <= 0 error
  # ------------------------------------------------------------
  if (!file.exists(snpfile)) {
    cat("  [WARN] SNP extract file not found, skipping:", snpfile, "\n")
    next
  }

  df_harm <- read.table(sumstatfile1, header = TRUE, stringsAsFactors = FALSE)
  snp_extract <- scan(snpfile, what = "character", quiet = TRUE)
  n_overlap <- sum(df_harm$SNP %in% snp_extract)

  cat("  Overlap between harmonized sumstats and extract SNPs:", n_overlap, "variants\n")
  if (n_overlap <= 0) {
    cat("  [ERROR] No SNP overlap with extract file for trait", ii,
        "MiXeR would fail (arg <= 0). Skipping this trait.\n")
    next
  }

  link <- paste0(
    "python3 ", dirPython, " fit1 ",
    "--trait1-file ", sumstatfile1, " ",
    "--out ", outfile, " ",
    "--extract ", snpfile, " ",
    "--bim-file ", bimfile, " ",
    "--ld-file ", ldfile, " ",
    "--lib ", lib
  )
  cat("  Running MiXeR command:\n", link, "\n")
  system(link)
}

cat("\nJob finished for sim", sim, "scenario", scenarioIndex, "\n")
