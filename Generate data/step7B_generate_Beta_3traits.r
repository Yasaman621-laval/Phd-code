#!/usr/bin/env Rscript
rm(list = ls())
library(foreach)
library(iterators)
library(doParallel)

# -------------------------------------------------------------------
# 1. Read array-task index
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0)
  stop("Usage: Rscript stage4_parallel_gwas.R <task_id>")

task_id <- as.numeric(args[1])
cat("Array task ID =", task_id, "\n")

# -------------------------------------------------------------------
# 2. Parallel setup
# -------------------------------------------------------------------
ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ncores) || ncores <= 0)
  ncores <- parallel::detectCores()

registerDoParallel(cores = ncores)
cat("Using", ncores, "cores (workers =", getDoParWorkers(), ")\n")

# -------------------------------------------------------------------
# 3. Fixed parameters (from your original Stage-4)
# -------------------------------------------------------------------
nsim   <- 10
trait_ids <- 1:3
pop    <- "EUR"
model  <- "linear"

dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype"
dirTraits <- file.path(dirBase, "three_traits")

refAllelefile <- file.path(dirBase, "refereceAlleleSNPs.txt")
usedsnpfile   <- file.path(dirBase, "SNPs_allChrs.txt")

dataname <- file.path(dirBase, paste0(pop, "AllChrs_bedformat"))

# -------------------------------------------------------------------
# 4. Identify scenario folders automatically (same pattern as your template)
# -------------------------------------------------------------------
scenario_dirs <- list.dirs(dirTraits, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0)
  stop("No scenario folder found under: ", dirTraits)

if (task_id < 1 || task_id > length(scenario_dirs))
  stop("task_id must be between 1 and ", length(scenario_dirs))

scenario_dir <- scenario_dirs[task_id]
scenario_name <- basename(scenario_dir)

cat("\n========================================================\n")
cat("Running Stage-4 GWAS for scenario:", scenario_name, "\n")
cat("Scenario path:", scenario_dir, "\n")
cat("========================================================\n\n")

# -------------------------------------------------------------------
# 5. Check phenotype files exist
# -------------------------------------------------------------------
for (tr in trait_ids) {
  test_pheno <- file.path(scenario_dir,
                          paste0("training.pheno_", pop, "_trait", tr, "_sim1"))

  if (!file.exists(test_pheno)) {
    stop("Missing phenotype files for trait ", tr,
         ". Run Stage-2 first.\nMissing file: ", test_pheno)
  }
}

# -------------------------------------------------------------------
# 6. Run GWAS inside scenario — PARALLELIZED over (trait × sim)
# -------------------------------------------------------------------
foreach(trait = trait_ids,
        .combine = "c") %:%        # nested foreach ? replaces your double-loop
foreach(sim = 1:nsim,
        .combine = "c") %dopar% {

  outname <- file.path(
    scenario_dir,
    paste0(pop, "_trait", trait, "_AllChrs_bedformat_sim", sim)
  )

  lfile <- paste0(outname, ".assoc.", model)
 # if (file.exists(lfile)) {
 #  cat("Skipping existing:", lfile, "\n")
   # return(NULL)
  #}

  phenofile <- file.path(
    scenario_dir,
    paste0("training.pheno_", pop, "_trait", trait, "_sim", sim)
  )

  if (!file.exists(phenofile)) {
    cat("Missing pheno file for trait", trait, "sim", sim, "\n")
    return(NULL)
  }

  plink_bin <- Sys.which("plink")
  if (!nzchar(plink_bin))
    stop("PLINK is not loaded. Load PLINK module first.")

  cmd <- paste(
    plink_bin,
    "--threads 1 --noweb --allow-no-sex",
    "--bfile", dataname,
    "--reference-allele", refAllelefile,
    "--extract", usedsnpfile,
    paste0("--", model),
    "--pheno", phenofile,
    "--out", outname
  )

  cat("Running:", cmd, "\n")
  status <- system(cmd)

  if (status != 0)
    stop("PLINK failed with status", status, "for", outname)

  cat("Completed:", outname, "\n")
  NULL
}

cat("\n? Stage-4 GWAS completed for scenario:", scenario_name, "\n")
