rm(list=ls())
library(MASS)
library(mvtnorm)
library(gtools)

#------------------------------------------------------------
# 1. Base setup
#------------------------------------------------------------
popvec <- c("AFR","EAS","EUR")

dirBase <- "/lustre10/scratch/yatah3/yatah3/simulation/SimuGenotype/"
#system(paste0("mkdir -p \"", dirBase, "\""))

totalsfile <- file.path(dirBase, "SNPs_allChrs.txt")
AllSNPs <- read.table(totalsfile, as.is=TRUE)[,1]
n.causal.snpPer <- 0.1
n.causal.snp <- ceiling(length(AllSNPs) * n.causal.snpPer)

#------------------------------------------------------------
# 2. Realistic simulation grid (EUR has largest training N)
#------------------------------------------------------------
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

#------------------------------------------------------------
# 3. Container for all results
#------------------------------------------------------------
summary_results <- data.frame()

#------------------------------------------------------------
# 4. Simulation loop
#------------------------------------------------------------
for (i in 1:length(rho_grid)) {

  rho <- rho_grid[i]
  hsq <- hsq_grid[i]
  TrainingNsam <- train_grid[[i]]
  valNsam <- val_grid[[i]]

  cat("\n==============================\n")
  cat("Scenario", i, ": hsq =", hsq, " rho =", rho,
      " Train =", paste(TrainingNsam, collapse = ","), "\n")
  cat("==============================\n")

  # Output directory for this setting
  subdir <- file.path(dirBase,
                      paste0("sim_hsq", hsq, "_rho", rho,
                             "_train", paste(TrainingNsam, collapse="-")))
  system(paste0("mkdir -p \"", subdir, "\""))

  #--------------------------------------------------------
  # Step 1. Generate causal SNPs and betas
  #--------------------------------------------------------
  set.seed(200 + i)
  causalSNPs <- sample(AllSNPs, n.causal.snp, replace=FALSE)

  var_beta <- hsq / n.causal.snp
  sigmause <- matrix(var_beta, 3, 3)
  for (ii in 1:3) {
    for (jj in 1:3) {
      if (ii != jj) sigmause[ii, jj] <- sigmause[ii, jj] * rho
    }
  }

  bvectors <- rmvnorm(n.causal.snp,
                      mean = rep(0, 3),
                      sigma = sigmause)

  tau <- quantile(abs(bvectors), 0.2)
  AllBetaMatrix <- matrix(0, nrow=length(AllSNPs), ncol=3)
  causal_positions <- match(causalSNPs, AllSNPs)
  AllBetaMatrix[causal_positions, ] <- bvectors

  vec1 <- AllBetaMatrix[, 1]
  vec2 <- AllBetaMatrix[, 2]
  vec3 <- AllBetaMatrix[, 3]

  vec1 <- abs(vec1) > tau
  vec2 <- abs(vec2) > tau
  vec3 <- abs(vec3) > tau
   
pi_0   <- mean(vec1 == 0 & vec2 == 0 & vec3 == 0)
pi_1   <- mean(vec1 != 0 & vec2 == 0 & vec3 == 0)
pi_2   <- mean(vec1 == 0 & vec2 != 0 & vec3 == 0)
pi_3   <- mean(vec1 == 0 & vec2 == 0 & vec3 != 0)
pi_12  <- mean(vec1 != 0 & vec2 != 0 & vec3 == 0)
pi_13  <- mean(vec1 != 0 & vec2 == 0 & vec3 != 0)
pi_23  <- mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
pi_123 <- mean(vec1 != 0 & vec2 != 0 & vec3 != 0)

pi_values <- c(pi_0, pi_1, pi_2, pi_3, pi_12, pi_13, pi_23, pi_123)
  print(pi_values)

  save(causalSNPs, bvectors, sigmause, var_beta, hsq,
       pi_values, file=file.path(subdir, "trueSNP_setting.RData"))

  #--------------------------------------------------------
  # Step 3. Simulate phenotype & compute PRS correlations
  #--------------------------------------------------------
  R2_vec <- numeric(3)
  names(R2_vec) <- popvec

  for (p in 1:3) {
    pop <- popvec[p]
    scorefile <- file.path(subdir, paste0("score_", pop))
    maffile <- file.path(dirBase, paste0(pop, "AllChrs_bedformat.frq"))
    dataname <- file.path(dirBase, paste0(pop, "AllChrs_bedformat"))
    tempfile <- file.path(subdir, paste0("temp_", pop))

    if (!file.exists(maffile)) {
      stop(paste("Missing MAF file for", pop, ":", maffile))
    }

    MAFdata0 <- read.table(maffile, as.is=TRUE, header=TRUE)
    x <- MAFdata0[match(causalSNPs, MAFdata0$SNP), c("SNP", "A1", "MAF")]
    score0 <- cbind(x[,1:2], bvectors[,p])
    write.table(score0, file=scorefile,
                row.names=FALSE, col.names=FALSE, quote=FALSE)

    # Run PLINK scoring
    cmd <- paste0("plink --bfile ", dataname,
                  " --noweb --allow-no-sex --score ",
                  scorefile, " --out ", tempfile)
    system(cmd)

    prof_path <- paste0(tempfile, ".profile")
    if (file.exists(prof_path)) {
      xb <- read.table(prof_path, as.is=TRUE, skip=1)[, c(1,2,6)]
      y <- xb
      y[,3] <- y[,3]/sqrt(var(y[,3])) * sqrt(sum(hsq)) +
               rnorm(nrow(y))*sqrt(1 - sum(hsq))
      R2_vec[p] <- cor(xb[,3], y[,3])^2
      cat(pop, ": R² =", round(R2_vec[p], 4), "\n")
    } else {
      warning(paste("PLINK output missing for", pop))
      R2_vec[p] <- NA
    }
  }

  #--------------------------------------------------------
  # Step 4. Append to summary
  #--------------------------------------------------------
  summary_results <- rbind(summary_results,
                           data.frame(
                             scenario = i,
                             rho = rho,
                             hsq = hsq,
                             Train_AFR = TrainingNsam[1],
                             Train_EAS = TrainingNsam[2],
                             Train_EUR = TrainingNsam[3],
                             R2_AFR = R2_vec["AFR"],
                             R2_EAS = R2_vec["EAS"],
                             R2_EUR = R2_vec["EUR"],
                             t(pi_values)
                           ))
}

#------------------------------------------------------------
# 5. Save summary
#------------------------------------------------------------
outfile <- file.path(dirBase, "SimulationSummary_realisticScenarios_realdata.csv")
write.csv(summary_results, outfile, row.names=FALSE)
cat("\n? Simulation complete. Summary saved to:\n", outfile, "\n")
