#!/usr/bin/env Rscript
rm(list = ls())

############################################################
# STEP 4 – PRS-CSx: Collect best R² for
#   - AFR (non-meta)
#   - EAS (non-meta)
#   - EUR–AFR (META)
#   - EUR–EAS (META)
# and plot them together (4 boxplots)
############################################################

########## Libraries ##########
library(dplyr)
library(ggplot2)

########## Read SLURM job input ##########
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing SLURM_ARRAY_TASK_ID")

jset <- as.numeric(args[1])

input_file <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/input_scenario.txt"
inputs     <- read.table(input_file, header = TRUE, as.is = TRUE)

scenarioIndex <- inputs[jset, 1]

cat(">>> scenarioIndex =", scenarioIndex, "\n")

########## General setup from your structure ##########
savename     <- "CSx-meta"
dirOutput    <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/"
path_to_python <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/PRScsx/"

popvec <- c("AFR", "EAS", "EUR")

# ============================================================
# 3. Scenario + TrainingNsam (same as your code)
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs)) stop("scenarioIndex exceeds available scenarios")

sc_dir  <- scenario_dirs[scenarioIndex]
dirSimuDataSet <- paste0(sc_dir, "/")
sc_name <- basename(sc_dir)

cat(">>> Using scenario:", sc_name, "\n")

# parse TrainingNsam from folder name suffix after "_train"
train_str       <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)
nums            <- as.numeric(
  sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x)))
)
TrainingNsam    <- nums

cat("    TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")

# ============================================================
# 4. phi values & general PRS-CSx setup
# ============================================================
phivec <- c(1e-6, 1e-4, 1e-2, 1)
nsim   <- 10

TrainIndex <- 3
popTrain   <- popvec[TrainIndex]  # EUR

# R2 output folders (non-meta and meta), consistent with your structure
R2output_nonmeta <- paste0(sc_dir, "/", savename, "/R2output/")
R2output_meta    <- paste0(sc_dir, "/", savename, "/R2output_meta/")

output_step4 <- paste0(sc_dir, "/output_step4_PRScsx/")
system(paste("mkdir -p", output_step4))

cat(">>> R2output_nonmeta:", R2output_nonmeta, "\n")
cat(">>> R2output_meta   :", R2output_meta, "\n")
cat(">>> output_step4    :", output_step4, "\n")

############################################################
# Utility 1: collect best R² across phi for NON-META (AFR/EAS)
#  - expects files:
#      R2output_nonmeta / (EUR_<pop>sim<sim>phiuse<phi>Test_results.RData)
############################################################
collect_R2_nonmeta <- function(pop, label_code) {

  cat(">>> Collecting NON-META PRS-CSx for pop =", pop, "\n")

  Start <- 1
  Allresult <- NULL

  for (sim in 1:nsim) {

    phi_list      <- c()
    R2_train_list <- c()   # EUR ? pop (validation)
    R2_test_list  <- c()   # pop ? pop (validation)
    nb_list       <- c()

    for (phiuse in phivec) {

      savecore <- paste0(popTrain, "_", pop,
                         "sim", sim, "phiuse", phiuse)

      train_file <- paste0(R2output_nonmeta, savecore, "Train_results.RData")
      test_file  <- paste0(R2output_nonmeta, savecore, "Test_results.RData")

      if (!file.exists(train_file) || !file.exists(test_file)) next

      #### A) EUR ? pop ---------------
      load(train_file)   # loads PreR2, numbetasvec
      eur_test_vec  <- as.numeric(PreR2[1, ])  # row 1 ? TEST
      eur_valid_vec <- as.numeric(PreR2[2, ])  # row 2 ? VALIDATION
      eur_nb_vec    <- as.numeric(numbetasvec)

      # SELECT f based on TEST (row 1)
      best_idx <- which.max(eur_test_vec)

      # FINAL EUR?pop R² from VALIDATION (row 2)
      eur_final_r2 <- eur_valid_vec[best_idx]
      eur_final_nb <- eur_nb_vec[best_idx]

      rm(PreR2, numbetasvec)

      #### B) pop ? pop ---------------
      load(test_file)
      pop_test_vec  <- as.numeric(PreR2[1, ])     # row 1 (not used)
      pop_valid_vec <- as.numeric(PreR2[2, ])     # row 2 ? final result
      pop_nb_vec    <- as.numeric(numbetasvec)

      pop_final_r2 <- pop_valid_vec[best_idx]     # apply SAME index
      pop_final_nb <- pop_nb_vec[best_idx]

      rm(PreR2, numbetasvec)

      #### Save model-level results
      phi_list      <- c(phi_list,      phiuse)
      R2_train_list <- c(R2_train_list, eur_final_r2)
      R2_test_list  <- c(R2_test_list,  pop_final_r2)
      nb_list       <- c(nb_list,       pop_final_nb)
    }

    if (length(phi_list) == 0) next

    #### Choose best f using EUR?pop VALIDATION R² (row 2)
    final_idx <- which.max(R2_train_list)

    temp <- c(
      phi_list[final_idx],         # selected f
      R2_train_list[final_idx],    # EUR?pop validation R²
      R2_test_list[final_idx],     # pop?pop validation R²
      nb_list[final_idx],          # # SNPs
      label_code
    )

    if (Start == 1) {
      Allresult <- temp
      Start <- 0
    } else {
      Allresult <- rbind(Allresult, temp)
    }
  }

  Allresult <- as.matrix(Allresult)
  colnames(Allresult) <- c(
    "phi_best",
    "EUR_to_POP_R2_valid",
    "POP_to_POP_R2_valid",
    "nBetas",
    "label_code"
  )

  return(Allresult)
}



############################################################
# Utility 2: collect best R² across phi for META (EUR–AFR/EAS)
#  - expects files:
#      R2output_meta / (EUR_<pop>sim<sim>phiuse<phi>_META_Meta_results.RData)
############################################################
collect_R2_meta <- function(pop, label_code) {

  cat(">>> Collecting META PRS-CSx for EUR–", pop, "\n")

  Start <- 1
  Allresult <- NULL

  for (sim in 1:nsim) {

    phi_list      <- c()
    r2_test_list  <- c()  # needed for phi-selection
    r2_val_list   <- c()
    nb_list       <- c()

    for (phiuse in phivec) {

      savecore_meta <- paste0(popTrain, "_", pop,
                              "sim", sim, "phiuse", phiuse, "_META")

      meta_file <- paste0(R2output_meta, savecore_meta, "_Meta_results.RData")

      if (!file.exists(meta_file)) next

      #### LOAD META RESULTS
      load(meta_file)   # loads PreR2 + numbetasvec

      test_vec  <- as.numeric(PreR2[1, ])   # row1 = test
      valid_vec <- as.numeric(PreR2[2, ])   # row2 = validation

      # choose f using test R²
      best_idx <- which.max(test_vec)

      # get selected test and validation R²
      best_test_r2  <- test_vec[best_idx]
      best_valid_r2 <- valid_vec[best_idx]
      best_nb       <- numbetasvec[best_idx]

      rm(PreR2, numbetasvec)

      #### Save results for each f
      phi_list      <- c(phi_list,      phiuse)
      r2_test_list  <- c(r2_test_list,  best_test_r2)
      r2_val_list   <- c(r2_val_list,   best_valid_r2)
      nb_list       <- c(nb_list,       best_nb)
    }

    # skip if empty sim
    if (length(phi_list) == 0) next

    #### select f using TEST R² (correct)
    best_idx_sim <- which.max(r2_test_list)

    temp <- c(
      phi_list[best_idx_sim],      # best phi
      r2_val_list[best_idx_sim],   # final R² = validation R²
      nb_list[best_idx_sim],       # # betas
      max(nb_list, na.rm = TRUE),  # max betas
      label_code
    )

    if (Start == 1) {
      Allresult <- temp
      Start <- 0
    } else {
      Allresult <- rbind(Allresult, temp)
    }
  }

  if (is.null(Allresult)) stop("No META results found.")

  Allresult <- as.matrix(Allresult)
  colnames(Allresult) <- c(
    "phi_best",
    "R2_valid",
    "nBetas_best",
    "nBetas_max",
    "label_code"
  )

  return(Allresult)
}


############################################################
# 1) Collect all result sets
############################################################

# NON-META
All_AFR_nonmeta <- collect_R2_nonmeta(pop = "AFR", label_code = 1)
All_EAS_nonmeta <- collect_R2_nonmeta(pop = "EAS", label_code = 2)

save(All_AFR_nonmeta,
     file = file.path(output_step4, "Allresult_PRScsx_AFR_nonmeta_full.RData"))
save(All_EAS_nonmeta,
     file = file.path(output_step4, "Allresult_PRScsx_EAS_nonmeta_full.RData"))

# META (EUR–AFR, EUR–EAS)
All_AFR_meta <- collect_R2_meta(pop = "AFR", label_code = 3)
All_EAS_meta <- collect_R2_meta(pop = "EAS", label_code = 4)

save(All_AFR_meta,
     file = file.path(output_step4, "Allresult_PRScsx_AFR_meta_full.RData"))
save(All_EAS_meta,
     file = file.path(output_step4, "Allresult_PRScsx_EAS_meta_full.RData"))

############################################################
# 2) Save 6 separate R2 vectors (for clarity)
############################################################

# NON-META: AFR
AFR_within_R2   <- All_AFR_nonmeta[, "POP_to_POP_R2_valid"]   # AFR?AFR
EUR_to_AFR_R2   <- All_AFR_nonmeta[, "EUR_to_POP_R2_valid"]   # EUR?AFR

save(AFR_within_R2,
     file = file.path(output_step4, "Allresult_PRScsx_AFR_within_nonmeta.RData"))
save(EUR_to_AFR_R2,
     file = file.path(output_step4, "Allresult_PRScsx_EUR_to_AFR_nonmeta.RData"))

# NON-META: EAS
EAS_within_R2   <- All_EAS_nonmeta[, "POP_to_POP_R2_valid"]   # EAS?EAS
EUR_to_EAS_R2   <- All_EAS_nonmeta[, "EUR_to_POP_R2_valid"]   # EUR?EAS

save(EAS_within_R2,
     file = file.path(output_step4, "Allresult_PRScsx_EAS_within_nonmeta.RData"))
save(EUR_to_EAS_R2,
     file = file.path(output_step4, "Allresult_PRScsx_EUR_to_EAS_nonmeta.RData"))

# META
AFR_meta_R2     <- All_AFR_meta[, "R2_valid"]
EAS_meta_R2     <- All_EAS_meta[, "R2_valid"]

save(AFR_meta_R2,
     file = file.path(output_step4, "Allresult_PRScsx_AFR_meta.RData"))
save(EAS_meta_R2,
     file = file.path(output_step4, "Allresult_PRScsx_EAS_meta.RData"))


