#!/usr/bin/env Rscript
# -------------------------------------------------------------------
# UTF-8 safe script - per simulation MiXeR pipeline
# -------------------------------------------------------------------
rm(list = ls())
library(dplyr)
library(stringr)

# -------------------------------------------------------------------
# 1. Parse job index (array ID)
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("? Missing argument: SLURM_ARRAY_TASK_ID")
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
dir1000G  <- paste0(dirBase, "/ld_mixer/")
output_sub_folder <- file.path(sc_dir, paste0("Mixer_sim", sim))
system(paste0("mkdir -p ", output_sub_folder))

# Mixer python + lib paths
dirPython <- "/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/precimed/mixer.py"
lib       <- "/home/yatah3/projects/def-thchlava/yatah3/Mixer/mixer/src/build/lib/libbgmg.so"
popvec <- c("AFR", "EAS", "EUR")

# Parse training sample sizes from scenario name
train_str <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)
nums <- as.numeric(sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x))))
TrainingNsam <- nums
cat(" TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")

# -------------------------------------------------------------------
# 5. Create summary statistics for this simulation only
# -------------------------------------------------------------------
GWASbetafile <- file.path(sc_dir, paste0("GWASbeta0Standard_allchrssim", sim, ".txt"))
if (!file.exists(GWASbetafile)) stop("? Missing GWAS file: ", GWASbetafile)
cat(" Reading:", GWASbetafile, "\n")

GWASbeta <- read.delim(GWASbetafile, header = TRUE, sep = "\t", quote = "",
                       fill = TRUE, strip.white = TRUE, comment.char = "",
                       stringsAsFactors = FALSE)

# Add A2 if missing
if (!"A2" %in% names(GWASbeta)) {
  GWASbeta <- GWASbeta %>%
    mutate(A2 = str_extract(SNP, "(?<=_[A-Za-z])_[A-Za-z]+$") %>% str_remove("^_")) %>%
    select(CHR, SNP, A1, A2,
           Zobs1, b1, SE1, p1,
           Zobs2, b2, SE2, p2,
           Zobs3, b3, SE3, p3)
}

# Split per population and write out per-sim files
for (ii in seq_along(popvec)) {
  cols_needed <- c("CHR", "SNP", "A1",
                   paste0("b", ii),
                   paste0("SE", ii),
                   paste0("p", ii))
  missing_cols <- setdiff(cols_needed, names(GWASbeta))
  if (length(missing_cols) > 0) {
    cat(" Missing columns for", popvec[ii], ":", paste(missing_cols, collapse = ", "), "\n")
    next
  }

  tmp <- GWASbeta[, cols_needed]
  names(tmp) <- c("CHR", "SNP", "A1", "BETA", "SE", "P")
  outfile <- file.path(sc_dir, paste0("summary_stats_pop", ii, "_sim", sim, ".txt"))
  write.table(tmp, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("  Created:", outfile, "\n")
}

# -------------------------------------------------------------------
# 6. Harmonize per-sim summary files
# -------------------------------------------------------------------
for (ii in 1:3) {
  infile  <- paste0(sc_dir, "/summary_stats_pop", ii, "_sim", sim, ".txt")
  outfile <- paste0(sc_dir, "/summary_stats_pop", ii, "_sim", sim, "_harmonized.txt")

  if (!file.exists(infile)) next
  df <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)
  names(df) <- toupper(names(df))

  df <- df %>%
    mutate(
      CHR = if ("CHR" %in% names(.)) CHR else as.numeric(sub(":.*", "", SNP)),
      BP  = as.numeric(sub(".*:(\\d+)_.*", "\\1", SNP)),
      A1  = toupper(A1),
      A2  = if ("A2" %in% names(.)) toupper(A2) else sub(".*_[A-Z]+_([A-Z]+)$", "\\1", SNP),
      Z   = ifelse(SE != 0, BETA / SE, 0),
      N   = TrainingNsam[ii]
    ) %>%
    select(SNP, CHR, BP, A1, A2, N, Z) %>%
    arrange(CHR, BP)

  write.table(df, outfile, quote = FALSE, sep = "\t", row.names = FALSE)
  cat("  Harmonized:", outfile, "\n")
}

# -------------------------------------------------------------------
# 7. Run MiXeR per population for this sim
# -------------------------------------------------------------------
for (ii in 1:3) {
  cat("\n===========================================================\n")
  cat(">>> Processing population:", popvec[ii], "for sim", sim, "\n")
  cat("===========================================================\n")

  sumstatfile1 <- paste0(sc_dir, "/summary_stats_pop", ii, "_sim", sim, "_harmonized.txt")
  outfile <- paste0(output_sub_folder, "/", popvec[ii], "_sim", sim, ".fit.rep")

  snpfile <- paste0(dir1000G, popvec[ii], "_allchr.prune_maf0p05_rand2M_r2p8.snps")
  bimfile <- paste0(sc_dir, "/Chr@/PlinkFormat/", popvec[ii], "chr@bedformat.bim")
  ldfile  <- paste0(dir1000G, "chr@.", popvec[ii], ".ld")

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


