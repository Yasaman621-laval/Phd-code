rm(list=ls())
library(MASS)
library(mvtnorm)
library(gtools)

#------------------------------------------------------------
# 1. Base setup
#------------------------------------------------------------
popvec <- "EUR"
trait_ids <- 1:3
nsim <- 10                     # Number of replications

dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"

# Dataset folders
dirSimuDataSet <- paste0(dirBase, "/three_traits/")


# SNP list
totalsfile <- paste0(dirBase, "SNPs_allChrs.txt")
AllSNPs <- read.table(totalsfile, as.is=TRUE)[,1]

#------------------------------------------------------------
# 2. Scenario grid (must match Stage-1)
#------------------------------------------------------------
rho_grid <- c(0.2, 0.4, 0.6, 0.8, 0.9)
hsq_grid <- c(0.3, 0.5, 0.6, 0.8, 0.9)

train_grid <- list(50000, 80000, 100000, 150000, 200000)
val_grid   <- list(5000,  10000, 15000, 20000, 25000)

#------------------------------------------------------------
# 3. Main simulation loop (same style as Stage-1)
#------------------------------------------------------------
for (i in 1:length(rho_grid)) {

  rho <- rho_grid[i]
  hsq <- hsq_grid[i]
  valNsam <- val_grid[[i]]
  TrainingNsam <- train_grid[[i]]

  cat("\n=====================================\n")
  cat("Stage-2 Scenario", i,
      "| hsq =", hsq,
      "| rho =", rho,
      "| Train_EUR =", TrainingNsam,
      "| Val_EUR =", valNsam, "\n")
  cat("=====================================\n")

  # Scenario-specific folder
  subdir <- paste0(
    dirSimuDataSet,
    "sim_hsq", hsq,
    "_rho", rho,
    "_trainEUR", TrainingNsam, "/"
  )

  trueSNPfile <- paste0(subdir, "trueSNP_setting.RData")
  load(trueSNPfile)       # Loads: causalSNPs, bvectors, sigmause, hsq, pi_values

  pop <- popvec
  dataname <- paste0(dirBase, pop, "AllChrs_bedformat")

  fam <- read.table(paste0(dataname, ".fam"), as.is=TRUE)
  samids <- 1:nrow(fam)

  #----------------------------------------------------------
  # 4. Repeat 20 phenotype simulations per scenario
  #----------------------------------------------------------
  for (sim in 1:nsim) {

    cat("\n  ---- Running sim", sim, "for scenario", i, "----\n")

    for (trait in trait_ids) {

      # Temporary PLINK output from Stage-1
      tempfile <- paste0(subdir, "temp_", pop, "_trait", trait, ".profile")

      if (!file.exists(tempfile)) {
        stop("Missing PLINK profile file: ", tempfile)
      }

      xb_y <- read.table(tempfile, as.is=TRUE, skip=1)[, c(1,2,6)]
      xb <- xb_y
      y  <- xb_y

      # phenotype simulation (same as Stage-1)
      y[,3] <- y[,3]/sqrt(var(y[,3])) * sqrt(sum(hsq)) +
               rnorm(nrow(y)) * sqrt(1 - sum(hsq))

      cat("    Trait", trait, "| R² =", round(cor(xb[,3], y[,3])^2, 4), "\n")

      #------------------------------------------------------
      # 5. Construct Test / Validation / Discovery splits
      #------------------------------------------------------

      test_ids  <- samids[1:(valNsam/2)]
      vali_ids  <- samids[((valNsam/2)+1):valNsam]
      disc_ids  <- samids[(valNsam+1):nrow(fam)]

      test_df <- fam[test_ids, 1:2]
      vali_df <- fam[vali_ids, 1:2]
      disc_df <- fam[disc_ids, 1:2]

      testvali_df <- fam[c(test_ids, vali_ids), 1:2]

      #------------------------------------------------------
      # 6. Write ID files
      #------------------------------------------------------
      testfile <- paste0(subdir,"TestDid", pop, "_trait", trait, "_sim", sim, ".txt")
      valfile  <- paste0(subdir,"ValiDid", pop, "_trait", trait, "_sim", sim, ".txt")
      discfile <- paste0(subdir,"Discoveryid", pop, "_trait", trait, "_sim", sim, ".txt")
      tvfile   <- paste0(subdir,"Test_ValiDid", pop, "_trait", trait, "_sim", sim, ".txt")

      write.table(test_df, file=testfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
      write.table(vali_df, file=valfile, quote=FALSE, row.names=FALSE, col.names=FALSE)
      write.table(disc_df,file=discfile,quote=FALSE, row.names=FALSE, col.names=FALSE)
      write.table(testvali_df,file=tvfile,quote=FALSE,row.names=FALSE,col.names=FALSE)

      #------------------------------------------------------
      # 7. Write phenotype files (train + test+val)
      #------------------------------------------------------
      ytestvalfile <- paste0(subdir, "test_validation.pheno_", pop,
                             "_trait", trait, "_sim", sim)
      ydiscfile <- paste0(subdir, "training.pheno_", pop,
                           "_trait", trait, "_sim", sim)

      write.table(y[c(test_ids, vali_ids),],
                  file=ytestvalfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

      write.table(y[disc_ids,],
                  file=ydiscfile, row.names=FALSE, col.names=FALSE, quote=FALSE)

    } # end trait loop
  } # end sim loop
} # end scenario loop

cat("\n? Stage-2 complete for all scenarios.\n")
