#!/usr/bin/env Rscript
rm(list = ls())


  library(SummaryLasso)
  library(dplyr)
  library(data.table)
  library(stringr)

# ============================================================
# 1. Parse job index (SLURM array ID)
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("? Missing argument: SLURM_ARRAY_TASK_ID")
jset <- as.numeric(args[1])

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project2/step1/"
inputs <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE, as.is = TRUE)

sim           <- inputs[jset, 1]
runIndex      <- inputs[jset, 2]
chrIndex      <- inputs[jset, 3]
scenarioIndex <- inputs[jset, 4]

cat(">>> sim =", sim, " runIndex =", runIndex, " chr =", chrIndex, " scenario =", scenarioIndex, "\n")

# ============================================================
# 2. Directories & scenario detection
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs)) stop("scenarioIndex exceeds available scenarios (", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)
cat(">>> Using scenario:", sc_name, "\n")

# ============================================================
# 3. Parse Training sample sizes
# ============================================================
train_str       <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)
TrainingNsam    <- as.numeric(sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x))))
cat(">>> TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")

# ============================================================
# 4. General settings
# ============================================================
savename     <- "new_densityU"
warmStart    <- 1
Zscale       <- 1
K            <- 2
singleStart  <- 1
numsnp       <- 400
rho_vec      <- c(seq(0,0.9,0.1),0.95)
groupnum     <- 2
nsim         <- 10
penalty = "RealmixLOG"
popvec       <- c("AFR","EAS","EUR")
functions_folder <- "/lustre09/project/6005709/yatah3/simulation/Rfunction/"

dirSimuData  <- paste0(sc_dir, "/")
output_sub_folder <- file.path(sc_dir, "Trans2-new-3/")


output_sub_folder_uDensity <- file.path(sc_dir, savename, "/")
system(paste0("mkdir -p ", output_sub_folder_uDensity))

# ============================================================
# 5. Load variance of effect sizes from new summary .csv files
# ============================================================
output_sub_folder0 <- paste0("/lustre09/project/6005709/yatah3/simulation/SimuGenotype/",sc_name,"/outputvisualonepop_sim",sim,"/")
var_beta <- numeric(length(popvec))

for (ii in 1:3) {
  finaloutputfit <- paste0(output_sub_folder0, popvec[ii], "_sim",sim,".csv")
  fit <- read.table(finaloutputfit, sep = ",", header = TRUE, as.is = TRUE)
  var_beta[ii] <- fit$sig2_beta[fit$Population == popvec[ii]]
}

# ============================================================
# 6. Choose pop-pair indices based on runIndex
# ============================================================
all_usedtrait <- c("1,3", "2,3")  # (AFR,EUR) or (EAS,EUR)
usedtrait     <- all_usedtrait[runIndex]
usedtraitsvec <- as.numeric(unlist(strsplit(usedtrait, ",")))
mainindex     <- usedtraitsvec[1]
iiIndexT      <- usedtraitsvec[2]
cat(">>> Trait pair:", paste(popvec[mainindex], popvec[iiIndexT], sep = "+"), "\n")

# ============================================================
# 7. Load Tau matrix for this pair
# ============================================================
Taufile <- paste0(sc_dir, "/", popvec[mainindex], "_", popvec[iiIndexT], "Tau_info_v2.RData")
if (!file.exists(Taufile)) stop("Taufile missing:", Taufile)
load(Taufile)  # loads AbsTauvec

# ============================================================
# 8. Load helper functions
# ============================================================
source(paste0(functions_folder, "AllPRS_Rfunctions.r"))
source(paste0(functions_folder, "PRS_utility.r"))
source(paste0(functions_folder, "Iterative_Rfunctions.r"))
source(paste0(functions_folder, "PlinkLD_transform.R"))

# ============================================================
# 9. Sigma2K vector (based on runIndex)
# ============================================================
if (runIndex == 1) {
  sigma2Kvec <- var_beta[c(1,3)]
} else {
  sigma2Kvec <- var_beta[c(2,3)]
}
cat(">>> sigma2Kvec =", paste(round(sigma2Kvec, 8), collapse = ", "), "\n")


 SharedPattern = permutations(n=2,r=K,v=c(0,1),repeats.allowed=T)
  num_alpha = nrow(SharedPattern)


# ========
#====================================================
# 10. Main loop over chromosomes and Tau grid
# ============================================================
for (chr in chrIndex) {
  print(paste0(">>> chr=", chr))
  
  for (tt in 1:ncol(AbsTauvec)) {
    print(paste0(">>> Tau column=", tt))
    tauuse <- AbsTauvec[sim, tt]
    
    subtau_saveoutfile <- paste0(output_sub_folder_uDensity, penalty, "_chr", chr, 
                                 "_usedtrait", usedtrait, "_warmStart", warmStart, 
                                 "_sim", sim, "_Zscale", Zscale, "_singleStart", singleStart, 
                                 "_tauuse", tauuse, "_DensityU.RData")
    
    dirSimuDataSetChr <- file.path(sc_dir, paste0("Chr", chr))
    dirPlinkFormat    <- file.path(dirSimuDataSetChr, "PlinkFormat")
    gwasfilename      <- paste0("GWASbetaStandard_allchrs", "sim", sim, ".txt")
    
    outnames_vec <- c()
    mafnames_vec <- c()
    for (ii in usedtraitsvec) {
      pop <- popvec[ii]
      base <- paste0(dirPlinkFormat, "/", pop, "chr", chr, "bedformat_numsnp", numsnp)
      outnames_vec <- c(outnames_vec, paste0(base, ".ld"))
    }
    
    Nvec <- TrainingNsam[usedtraitsvec]
    GWASbetafile <- file.path(dirSimuData, gwasfilename)
    
    summaryZoutput <- Get_summaryZ(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr)
    summaryZ       <- summaryZoutput[[1]]
    Allkeepsnps    <- summaryZoutput[[2]]
    rm(summaryZoutput)
    
    plinkLD <- Gen_All_lddata_simpleTrans(outnames_vec, Allkeepsnps, rcut = 0.3)
    
    # Generate Beta matrices
    savefile1 <- paste0(output_sub_folder, penalty,"chr", chr, "usedtrait", usedtrait, "warmStart", warmStart,
                        "sim", sim, "Zscale", Zscale, "tauuse", tauuse, "singleStart", singleStart, ".RData")
    
    BetaMatrixS <- Gen_One_BetaMatrix(savefile1, 2, 1)[[1]]
    BetaMatrix3 <- Gen_One_BetaMatrix(savefile1, 2, 2)[[1]]
    
    mat <- match(rownames(summaryZ), colnames(BetaMatrixS))
    BetaMatrixS <- BetaMatrixS[, mat]
    mat <- match(rownames(summaryZ), colnames(BetaMatrix3))
    BetaMatrix3 <- BetaMatrix3[, mat]
    
    # Combine Beta matrices
    AllBetaMatrix <- do.call(rbind, lapply(1:nrow(BetaMatrixS), function(i) c(BetaMatrixS[i,], BetaMatrix3[i,])))
    
    # Generate es_alphaMatrix
    es_alphaMatrix <- matrix(nrow = nrow(AllBetaMatrix), ncol = 4)
    for (i in 1:nrow(AllBetaMatrix)) {
      non_zero <- AllBetaMatrix[i, ][AllBetaMatrix[i, ] != 0]
      if (length(non_zero) == 0) {
        result2 <- c(1, 0, 0, 0)
      } else {
        half_non_zero <- floor(length(non_zero) / length(Nvec))
        reshaped_beta <- matrix(0, nrow = ncol(BetaMatrixS), ncol = length(Nvec))
        reshaped_beta[1:half_non_zero, 1] <- non_zero[1:half_non_zero]
        reshaped_beta[1:(length(non_zero) - half_non_zero), length(Nvec)] <- non_zero[(half_non_zero + 1):length(non_zero)]
        vec1 <- reshaped_beta[, 1]; vec2 <- reshaped_beta[, 2]
        result2 <- c(mean(vec1 == 0 & vec2 == 0),
                     mean(vec1 == 0 & vec2 != 0),
                     mean(vec1 != 0 & vec2 == 0),
                     mean(vec1 != 0 & vec2 != 0))
      }
      es_alphaMatrix[i, ] <- result2
    }
    
    rm(BetaMatrix3, BetaMatrixS)
    gc()
    
    # Posterior computation for all rho values
    Allrho_WBMatrix_list <- list()
    for (rr in seq_along(rho_vec)) {
      rho <- rho_vec[rr]
      rhoMat <- matrix(rho, K, K)
      sigma2K_allAlpha_List <- Create_sigma2K_allAlpha_List(K, sigma2Kvec, rhoMat)

    output <- transConditionalU(summaryZ = summaryZ,Nvec = Nvec,JointBmatrix = AllBetaMatrix, Zcov = diag(K),
  SDvec = NULL,
  Zscale = 1,
  sigma2K_allAlpha_List = sigma2K_allAlpha_List,
  es_alphaMatrix = es_alphaMatrix,
  plinkLD = plinkLD
)
     
      
      SNPnames <- rownames(summaryZ)
      P <- length(SNPnames)
      inputs1 <- cbind(output[[1]], output[[4]])
      wBMatrix1 <- apply(inputs1, 1, gen_sPostB, P, num_alpha)
      rownames(wBMatrix1) <- SNPnames
      inputs2 <- cbind(output[[1]], output[[5]])
      wBMatrix2 <- apply(inputs2, 1, gen_sPostB, P, num_alpha)
      rownames(wBMatrix2) <- SNPnames
      
      wBoutput <- list(wBMatrix1, wBMatrix2)
      Allrho_WBMatrix_list[[rr]] <- wBoutput
      rm(wBoutput, wBMatrix1, wBMatrix2, output)
      gc()
    }
    
    save(Allrho_WBMatrix_list, file = subtau_saveoutfile)
    rm(Allrho_WBMatrix_list, es_alphaMatrix, AllBetaMatrix)
    gc()
  }
}
cat("Completed sim", sim, "runIndex", runIndex, "chr", chrIndex, "scenario", scenarioIndex, "\n")
