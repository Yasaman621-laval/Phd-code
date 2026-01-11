rm(list=ls())
library(MASS)
library(mvtnorm)
library(gtools)

#------------------------------------------------------------
# 1. Base setup
#------------------------------------------------------------
popvec <- "EUR"            # Only EUR population
trait_ids <- 1:3           # Three traits

dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"

dirSimuDataSet <- file.path(dirBase, "three_traits")

system(paste0("mkdir -p \"", dirSimuRoot, "\""))
system(paste0("mkdir -p \"", dirSimuDataSet, "\""))

totalsfile <- file.path(dirBase, "SNPs_allChrs.txt")
AllSNPs <- read.table(totalsfile, as.is = TRUE)[,1]

n.causal.snpPer <- 0.1
n.causal.snp <- ceiling(length(AllSNPs) * n.causal.snpPer)
cat("Total SNPs:", length(AllSNPs), "  Causal:", n.causal.snp, "\n")

#------------------------------------------------------------
# 2. Realistic simulation grid (5 scenarios)
#------------------------------------------------------------
rho_grid <- c(0.2, 0.4, 0.6, 0.8, 0.9)
hsq_grid <- c(0.3, 0.5, 0.6, 0.8, 0.9)

train_grid <- list(
  50000,   # Scenario A
  80000,   # Scenario B
  100000,  # Scenario C
  150000,  # Scenario D
  200000   # Scenario E
)

val_grid <- list(
  5000,
  10000,
  15000,
  20000,
  25000
)

#------------------------------------------------------------
# 3. Container for summary results
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

  cat("\n=====================================\n")
  cat("Scenario", i,
      "| hsq =", hsq,
      "| rho =", rho,
      "| Train_EUR =", TrainingNsam,
      "| Val_EUR =", valNsam, "\n")
  cat("=====================================\n")

  # Scenario directory
  subdir <- file.path(
    dirSimuDataSet,
    paste0("sim_hsq", hsq,
           "_rho", rho,
           "_trainEUR", TrainingNsam)
  )
  system(paste0("mkdir -p \"", subdir, "\""))

  #--------------------------------------------------------
  # Step 1: generate causal SNPs and effect sizes
  #--------------------------------------------------------
  set.seed(200 + i)
  causalSNPs <- sample(AllSNPs, n.causal.snp, replace = FALSE)

  var_beta <- hsq / n.causal.snp
  sigmause <- matrix(var_beta, 3, 3)

  for (ii in 1:3) {
    for (jj in 1:3) {
      if (ii != jj) sigmause[ii, jj] <- sigmause[ii, jj] * rho
    }
  }

  cat("Sigma matrix:\n")
  print(sigmause)

  bvectors <- rmvnorm(
    n = n.causal.snp,
    mean = rep(0, 3),
    sigma = sigmause
  )

  # Threshold for causal classification
  tau <- quantile(abs(bvectors), 0.2)

  AllBetaMatrix <- matrix(0, nrow = length(AllSNPs), ncol = 3)
  causal_positions <- match(causalSNPs, AllSNPs)
  AllBetaMatrix[causal_positions, ] <- bvectors

  vec1 <- AllBetaMatrix[,1]
  vec2 <- AllBetaMatrix[,2]
  vec3 <- AllBetaMatrix[,3]

  causal_in1 <- abs(vec1) > tau
  causal_in2 <- abs(vec2) > tau
  causal_in3 <- abs(vec3) > tau

  pi_values <- c(
    pi_0   = mean(!causal_in1 & !causal_in2 & !causal_in3),
    pi_1   = mean(causal_in1 & !causal_in2 & !causal_in3),
    pi_2   = mean(!causal_in1 & causal_in2 & !causal_in3),
    pi_3   = mean(!causal_in1 & !causal_in2 & causal_in3),
    pi_12  = mean(causal_in1 & causal_in2 & !causal_in3),
    pi_13  = mean(causal_in1 & !causal_in2 & causal_in3),
    pi_23  = mean(!causal_in1 & causal_in2 & causal_in3),
    pi_123 = mean(causal_in1 & causal_in2 & causal_in3)
  )

  print(pi_values)

  trueSNPfile <- file.path(subdir, "trueSNP_setting.RData")
  save(causalSNPs, bvectors, sigmause, var_beta, hsq, pi_values,
       file = trueSNPfile)

  load(trueSNPfile)

  #--------------------------------------------------------
  # Step 2: simulate phenotypes + compute R² for each trait
  #--------------------------------------------------------
  R2_vec <- numeric(length(trait_ids))
  names(R2_vec) <- paste0("trait", trait_ids)

  pop <- popvec

  for (trait in trait_ids) {

    scorefile <- file.path(subdir, paste0("score_", pop, "_trait", trait))

    if (!file.exists(scorefile)) {

      maffile <- file.path(dirBase,
                           paste0(pop, "AllChrs_bedformat.frq"))
      if (!file.exists(maffile)) stop("Missing MAF file: ", maffile)

      MAFdata0 <- read.table(maffile, as.is = TRUE, header = TRUE)
      x <- MAFdata0[match(causalSNPs, MAFdata0$SNP), c("SNP","A1","MAF")]

      dataname <- file.path(dirBase,
                            paste0(pop, "AllChrs_bedformat"))

      score0 <- cbind(x[,1:2], bvectors[, trait])
      write.table(score0, file = scorefile,
                  row.names = FALSE, col.names = FALSE, quote = FALSE)

      tempfile <- file.path(subdir,
                            paste0("temp_", pop, "_trait", trait))

      cmd <- paste0(
        "plink --bfile ", dataname,
        " --noweb --allow-no-sex",
        " --score ", scorefile,
        " --out ", tempfile
      )
      system(cmd)

      prof_path <- paste0(tempfile, ".profile")
      if (!file.exists(prof_path)) {
        R2_vec[trait] <- NA
        next
      }

      xb_y <- read.table(prof_path, as.is = TRUE, skip = 1)[, c(1,2,6)]
      y <- xb_y
      y[,3] <- y[,3] / sqrt(var(y[,3])) *
               sqrt(sum(hsq)) +
               rnorm(nrow(y)) * sqrt(1 - sum(hsq))

      R2_vec[trait] <- cor(xb_y[,3], y[,3])^2
      cat(pop, "Trait", trait, ": R² =", round(R2_vec[trait], 4), "\n")
    }
  }

  #--------------------------------------------------------
  # Step 3: append scenario summary
  #--------------------------------------------------------
  summary_results <- rbind(
    summary_results,
    data.frame(
      scenario = i,
      rho = rho,
      hsq = hsq,
      Train_EUR = TrainingNsam,
      Val_EUR = valNsam,
      R2_trait1 = R2_vec["trait1"],
      R2_trait2 = R2_vec["trait2"],
      R2_trait3 = R2_vec["trait3"],
      t(pi_values)
    )
  )
}

#------------------------------------------------------------
# 5. Save summary
#------------------------------------------------------------
outfile <- file.path(dirSimuDataSet, "SimulationSummary_threeTraits.csv")
write.csv(summary_results, outfile, row.names = FALSE)

cat("\n? Simulation complete. Summary saved to:\n", outfile, "\n")
