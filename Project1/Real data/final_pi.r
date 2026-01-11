# =============================================================
# Compute global es_alphaMatrix across all chromosomes
# Using the same structure as the single-chr computation
# =============================================================
rm(list=ls())
library(dplyr)

# -------------------------------------------------------------
# PARAMETERS
# -------------------------------------------------------------
penalty <- "RealmixLOG"
savename <- "new_densityU_3pops_project1"
warmStart <- 1
Zscale <- 1
singleStart <- 1
K <- 3

dirOutput <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/project1/"
output_sub_folder <- "/lustre06/project/6005709/yatah3/real-data/project1/Trans2-new-3-out/"
output_sub_folder_uDensity <- paste0("/lustre06/project/6005709/yatah3/real-data/project1/", savename, "/")

popvec <- c("AFR","EAS","EUR")

# -------------------------------------------------------------
# LOAD Tau file and get target tau, rho, etc.
# -------------------------------------------------------------
Taufile <- paste0(dirOutput, popvec[1], "_", popvec[2], "-", popvec[3], "Tau_info_v2.RData")
load(Taufile)
if (!exists("AbsTauvec")) stop("AbsTauvec not found in Tau file")
tauuse <- AbsTauvec[which.min(abs(AbsTauvec - 7.5e-05))]
rho_vec <- c(seq(0,0.9,0.1),0.95)
rho_target <- 0.95

cat("Using rho =", rho_target, "tau =", tauuse, "\n")

  
functionsfolder  <- "/lustre06/project/6005709/yatah3/real-data/Rfunction/"

  source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))
  source(paste0(functionsfolder,"PRS_utility.r"))
  source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
  source(paste0(functionsfolder,"PlinkLD_transform.R"))





# -------------------------------------------------------------
# Combine BetaMatrixS / BetaMatrix3 / BetaMatrixE across all chr
# -------------------------------------------------------------

# -------------------------------------------------------------
# Initialize combined matrices
# -------------------------------------------------------------
BetaMatrixS_all <- NULL
BetaMatrix3_all <- NULL
BetaMatrixE_all <- NULL

# -------------------------------------------------------------
# Combine Beta matrices across chromosomes
# -------------------------------------------------------------
for (chr in 1:22) {
  cat("Reading Beta matrices for chr", chr, "...\n")

  savefile1 <- paste0(
    output_sub_folder, penalty,
    "chr", chr, "_3traits",
    "warmStart", warmStart,
    "Zscale", Zscale,
    "tauuse", tauuse,
    "singleStart", singleStart, ".RData"
  )

  if (!file.exists(savefile1)) {
    warning("Missing file for chr=", chr)
    next
  }

  # Extract the three population beta matrices
  BetaMatrixS <- Gen_One_BetaMatrix(savefile1, 3, 1)[[1]]
  BetaMatrix3 <- Gen_One_BetaMatrix(savefile1, 3, 2)[[1]]
  BetaMatrixE <- Gen_One_BetaMatrix(savefile1, 3, 3)[[1]]

  # Combine column-wise across chromosomes
  if (chr == 1) {
    BetaMatrixS_all <- BetaMatrixS
    BetaMatrix3_all <- BetaMatrix3
    BetaMatrixE_all <- BetaMatrixE
  } else {
    BetaMatrixS_all <- cbind(BetaMatrixS_all, BetaMatrixS)
    BetaMatrix3_all <- cbind(BetaMatrix3_all, BetaMatrix3)
    BetaMatrixE_all <- cbind(BetaMatrixE_all, BetaMatrixE)
  }
}

cat("Final combined dimensions:\n")
cat("S:", nrow(BetaMatrixS_all), "x", ncol(BetaMatrixS_all), "\n")
cat("3:", nrow(BetaMatrix3_all), "x", ncol(BetaMatrix3_all), "\n")
cat("E:", nrow(BetaMatrixE_all), "x", ncol(BetaMatrixE_all), "\n")


# -------------------------------------------------------------
# Compute global es_alphaMatrix exactly like original structure
# -------------------------------------------------------------
es_alphaMatrix <- matrix(0, nrow = 100, ncol = 8)
colnames(es_alphaMatrix) <- c("pi_0","pi_1","pi_2","pi_3","pi_12","pi_13","pi_23","pi_123")

for (i in 1:100) {
  # Extract global vectors (row i from combined matrices)
  vec1 <- BetaMatrixS_all[i, ]
  vec2 <- BetaMatrix3_all[i, ]
  vec3 <- BetaMatrixE_all[i, ]
  
  # Calculate proportions (following your exact formula)
  pi_0  <- mean(vec1 == 0 & vec2 == 0 & vec3 == 0)
  pi_1  <- mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
  pi_2  <- mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
  pi_3  <- mean(vec1 != 0 & vec2 != 0 & vec3 == 0)
  pi_12 <- mean(vec1 != 0 & vec2 != 0 & vec3 == 0)
  pi_13 <- mean(vec1 != 0 & vec2 == 0 & vec3 != 0)
  pi_23 <- mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
  pi_123 <- 1 - (pi_0 + pi_1 + pi_2 + pi_3 + pi_12 + pi_13 + pi_23)
  
  es_alphaMatrix[i, ] <- c(pi_0, pi_1, pi_2, pi_3, pi_12, pi_13, pi_23, pi_123)
}

cat("\nGlobal es_alphaMatrix computed successfully.\n")
print(head(es_alphaMatrix))

# -------------------------------------------------------------
# Extract row 68 for target rho and tau
# -------------------------------------------------------------
row68 <- es_alphaMatrix[68, ]
row68_df <- as.data.frame(t(row68))
row68_df$rho <- rho_target
row68_df$tau <- tauuse

cat("\nRow 68 summary (global, across all chr):\n")
print(row68_df)

# -------------------------------------------------------------
# Save outputs
# -------------------------------------------------------------
write.table(
  es_alphaMatrix,
  file = paste0(dirOutput, "es_alphaMatrix_combined_allchr.txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  row68_df,
  file = paste0(dirOutput, "es_alphaMatrix_row68_rho", rho_target, "_tau", tauuse, ".txt"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("\nAll files saved in:", dirOutput, "\n")


#_________________________
#results
#_________________________

Row 68 summary (global, across all chr):
       pi_0      pi_1      pi_2      pi_3     pi_12    pi_13     pi_23
 0.2620375 0.0848138 0.0848138 0.1315958 0.1315958 0.122572 0.0848138
      pi_123  
      0.09775751
      
      rho     tau
  0.95 7.5e-05
