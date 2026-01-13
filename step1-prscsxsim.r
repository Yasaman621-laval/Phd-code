#!/usr/bin/env Rscript
rm(list = ls())

cat("\n===== PRS-CSx META RUN STARTED =====\n")

# ============================================================
# 1. Parse unified input format (same as step1 & step2)
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing SLURM_ARRAY_TASK_ID")

jset <- as.numeric(args[1])

input_file <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/input.txt"

inputs <- read.table(input_file, header = TRUE, as.is = TRUE)

sim           <- inputs[jset, 1]
runIndex      <- inputs[jset, 2]
phiIndex      <- inputs[jset, 3]   # maps to ttIndex style
scenarioIndex <- inputs[jset, 4]

cat(">>> sim =", sim, 
    " runIndex =", runIndex, 
    " phiIndex =", phiIndex, 
    " scenario =", scenarioIndex, "\n")

# ============================================================
# 2. Fixed settings
# ============================================================

savename <- "CSx-meta"

dirOutput <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/"

path_to_python <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/PRScsx/"


popvec <- c("AFR", "EAS", "EUR")

# ============================================================
# 3. Load summary stats for sample sizes
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs)) stop("scenarioIndex exceeds available scenarios (", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]

output_sub_folder <- paste0(sc_dir,"/", savename, "/")
system(paste0("mkdir -p ", output_sub_folder))


sc_name <- basename(sc_dir)
cat(">>> Using scenario:", sc_name, "\n")

# parse TrainingNsam from folder name suffix after "_train"
train_str       <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)
nums            <- as.numeric(sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x))))
TrainingNsam <- as.integer(nums)


cat("    TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")
# ============================================================
# 4. phi values (ttIndex maps here)
# ============================================================

phivec <- c(1e-6, 1e-4, 1e-2, 1)

if (phiIndex > length(phivec))
  stop("phiIndex exceeds available phi values")

phiuse <- phivec[phiIndex]
cat(">>> using phi =", phiuse, "\n")

# ============================================================
# 5. Determine which population pair (runIndex)
# ============================================================

TrainIndex <- 3   # EUR always training
popTrain   <- popvec[TrainIndex]

if (runIndex == 1) {
  TargetIndex <- 1   # AFR
} else if (runIndex == 2) {
  TargetIndex <- 2   # EAS
} else {
  stop("Invalid runIndex (must be 1 or 2)")
}

popT <- popvec[TargetIndex]

cat(">>> Running pair: ", popTrain, "+", popT, "\n")

# ============================================================
# 6. Check BIM alignment
# ============================================================
bim_prefix <- paste0(dirBase, "/", "EURAllChrs_bedformat")

if (!file.exists(paste0(bim_prefix, ".bim")))
  stop("Missing BIM: ", paste0(bim_prefix, ".bim"))

# ============================================================
# 7. Summary statistics for training & target
# ============================================================

sumstatTrainfile <- paste0(sc_dir, "/GWAS_", popTrain,"_",sim, ".tsv")
sumstatfile      <- paste0(sc_dir, "/GWAS_", popT,"_", sim, ".tsv")

# ============================================================
# 8. Build PRS-CSx meta command
# ============================================================

# 8. Build PRS-CSx meta command (clean version)

outname <- paste0("PRSCSx_", popTrain, "_", popT, "_phi", phiuse, "_sim", sim)
python_bin <- "~/scratch/ENmixer311/bin/python"

cmd <- paste(
  python_bin, 
  file.path(path_to_python, "PRScsx.py"),
  paste0("--ref_dir=", dirOutput),
  paste0("--bim_prefix=", bim_prefix),
  paste0("--sst_file=", sumstatTrainfile, ",", sumstatfile),
  paste0("--n_gwas=", TrainingNsam[TrainIndex], ",", TrainingNsam[TargetIndex]),
  paste0("--pop=", popTrain, ",", popT),
  "--chrom=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22",
  paste0("--phi=", phiuse),
  paste0("--out_dir=", output_sub_folder),
  paste0("--out_name=", outname),
  "--meta=TRUE"
)

system(cmd, wait = TRUE, intern = FALSE)



cat("\n===== PRS-CSx META RUN COMPLETE =====\n")
