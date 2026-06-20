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

dirBase <- "/lustre10/scratch/yatah3/yatah3/simulation/SimuGenotype/"
scenario_dirs_all <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs_all <- scenario_dirs_all[
  grepl("^sim_hsq", basename(scenario_dirs_all))
]

# Optional but important: sort to make scenario order stable
scenario_dirs_all <- sort(scenario_dirs_all)

cat("\nAll detected scenarios:\n")
print(data.frame(
  original_index = seq_along(scenario_dirs_all),
  scenario = basename(scenario_dirs_all)
))

# Keep only scenarios 1, 3, and 4
selected_scenarios <- c(1, 3, 4)

scenario_dirs <- scenario_dirs_all[selected_scenarios]

cat("\nSelected scenarios:\n")
print(data.frame(
  original_index = selected_scenarios,
  local_index = seq_along(scenario_dirs),
  scenario = basename(scenario_dirs)
))

if (length(scenario_dirs) == 0) {
  stop("No selected scenario folders found under: ", dirBase)
}

# -----------------------------------------------------------
# Scenario grid (must match the simulation step)
# -----------------------------------------------------------
rho_grid <- c(0.65, 0.80, 0.95)

hsq_grid <- c(0.20, 0.30, 0.40)

train_grid <- list(
  c(20000, 5000, 70000),
  c(20000, 5000, 70000),
  c(20000, 5000, 70000)
)

val_grid <- list(
  c(5000, 1000, 3000),
  c(5000, 1000, 3000),
  c(5000, 1000, 3000)
)
# -----------------------------------------------------------
# helper to find scenario index from folder name (based on rho & hsq)
# -----------------------------------------------------------
get_scenario_index <- function(folder_name, scenario_dirs) {
  match(folder_name, basename(scenario_dirs))
}

# -----------------------------------------------------------
# main constants (same as your old script)
# -----------------------------------------------------------
nsim <- 10

# -----------------------------------------------------------
# MAIN loop over scenario folders
# -----------------------------------------------------------
for (local_idx in seq_along(scenario_dirs)) {

  sc_dir <- scenario_dirs[local_idx]
  sc_name <- basename(sc_dir)

  original_idx <- selected_scenarios[local_idx]

  valNsam_vec <- val_grid[[local_idx]]
names(valNsam_vec) <- popvec

cat("\n==============================\n")
cat("Processing original Scenario", original_idx, ":", sc_name, "\n")
cat("Local grid index:", local_idx, "\n")
cat("Validation sizes (per pop): ",
    paste0(popvec, "=", valNsam_vec, collapse = ", "), "\n")
cat("==============================\n")

  # ---------------------------------------------------------
  # Load true SNP info (contains hsq) - scenario-specific
  # ---------------------------------------------------------
  trueSNPfile <- paste0(sc_dir, "/trueSNP_setting.RData")
  if (!file.exists(trueSNPfile)) stop("Missing: ", trueSNPfile)

  load(file = trueSNPfile)  # must load hsq

   for (iiIndex in c(1:3)) {

    pop <- popvec[iiIndex]

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

# ---- NEW (ADDED: separate test & validation) ----
ytestfile <- paste0(sc_dir, "/test.pheno_", pop, "sim", sim)
yvalfile  <- paste0(sc_dir, "/validation.pheno_", pop, "sim", sim)

# ---- WRITE FILES ----

# training (Discovery)
write.table(
  y[didp, ],
  file = ydisfile,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# combined test + validation (keep for compatibility)
write.table(
  y[c(tidp, vidp), ],
  file = ytestvalfile,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# -------------------------
# NEW: validation only
# -------------------------
write.table(
  y[vidp, ],
  file = yvalfile,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

# -------------------------
# NEW: test only
# -------------------------
write.table(
  y[tidp, ],
  file = ytestfile,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)
    }
  }
}

cat("\nDone: OLD-style deterministic splits + OLD filenames, but scenario-aware folders.\n")
