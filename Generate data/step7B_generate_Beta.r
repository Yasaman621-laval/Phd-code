#!/usr/bin/env Rscript
options(encoding = "UTF-8")
rm(list = ls())

library(foreach)
library(iterators)
library(doParallel)

# -------------------------------------------------------------------
# 1. Parse job array index argument
# -------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No argument supplied. Usage: Rscript step7B_generate_Beta_args.R <task_id>")
}

task_id <- as.numeric(args[1])
cat("Array task ID received:", task_id, "\n")

# -------------------------------------------------------------------
# 2. Parallel setup
# -------------------------------------------------------------------
ncores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
if (is.na(ncores) || ncores <= 0) {
  ncores <- parallel::detectCores(logical = FALSE)
}
registerDoParallel(cores = ncores)
cat("Using", ncores, "cores for parallel PLINK runs\n")

# -------------------------------------------------------------------
# 3. Fixed parameters
# -------------------------------------------------------------------
nsim <- 10
model <- "linear"
popvec <- c("AFR", "EAS", "EUR")

dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype"

scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) {
  stop(paste("No scenario folders found under:", dirBase))
}

if (task_id < 1 || task_id > length(scenario_dirs)) {
  stop(paste("Invalid task_id. It must be between 1 and", length(scenario_dirs)))
}

# Pick scenario
sc_dir <- scenario_dirs[task_id]
sc_name <- basename(sc_dir)
cat("Running GWAS for scenario:", sc_name, "\n")

# Reference and SNP files
refAllelefile <- file.path(dirBase, "refereceAlleleSNPs.txt")
usedsnpfile   <- file.path(dirBase, "SNPs_allChrs.txt")

if (!file.exists(refAllelefile)) {
  stop(paste("Missing refereceAlleleSNPs.txt in", dirBase))
}
if (!file.exists(usedsnpfile)) {
  stop(paste("Missing SNPs_allChrs.txt in", dirBase))
}

# -------------------------------------------------------------------
# 4. Run PLINK GWAS for each population x simulation replicate
# -------------------------------------------------------------------
GeneBetaIndex <- 1

if (GeneBetaIndex == 1) {

  for (iiIndex in seq_along(popvec)) {

    pop <- popvec[iiIndex]
    dataname <- file.path(dirBase, paste0(pop, "AllChrs_bedformat"))

    if (!file.exists(paste0(dataname, ".bim"))) {
      message("Missing genotype file for ", pop, " - skipping.")
      next
    }

    foreach(sim = 1:nsim, .combine = rbind) %dopar% {

      phenofile <- file.path(
        sc_dir,
        paste0("training.pheno_", pop, "_sim", sim)
      )

      if (!file.exists(phenofile)) {
        message("Missing phenotype file: ", phenofile)
        return(NULL)
      }

      outname <- file.path(
        sc_dir,
        paste0(pop, "_AllChrs_bedformat_sim", sim)
      )

      plink_bin <- Sys.which("plink")
      if (!nzchar(plink_bin)) {
        stop("PLINK not found in PATH")
      }

      lineR <- paste0(
        plink_bin,
        " --noweb --allow-no-sex",
        " --bfile ", dataname,
        " --reference-allele ", refAllelefile,
        " --extract ", usedsnpfile,
        " --", model,
        " --pheno ", phenofile,
        " --out ", outname
      )

      system(lineR)
      NULL
    }
  }
}

cat("Completed GWAS for scenario:", sc_name, "\n")
