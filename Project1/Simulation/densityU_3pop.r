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

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project1/step1/"
inputs <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE, as.is = TRUE)

sim           <- inputs[jset, 1]
chrIndex           <- inputs[jset, 2]
scenarioIndex <- inputs[jset, 3]

cat(">>> sim =", sim, " chr =", chrIndex, " scenario =", scenarioIndex, "\n")

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
savename     <- "new_densityU_3pop"
warmStart    <- 1
Zscale       <- 1
K            <- 3
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
# 7. Load Tau matrix for this pair
# ============================================================
Taufile <- paste0(sc_dir,"/",popvec[1],"_",popvec[2],"_",popvec[3],"Tau_info_v2.RData")
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

sigma2Kvec <- var_beta

cat(">>> sigma2Kvec =", paste(round(sigma2Kvec, 8), collapse = ", "), "\n")


SharedPattern <- as.matrix(
  expand.grid(rep(list(c(0,1)), K))
)

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
                                 "_3pop_", "_warmStart", warmStart, 
                                 "_sim", sim, "_Zscale", Zscale, "_singleStart", singleStart, 
                                 "_tauuse", tauuse, "_DensityU.RData")
    
    dirSimuDataSetChr <- file.path(sc_dir, paste0("Chr", chr))
    dirPlinkFormat    <- file.path(dirSimuDataSetChr, "PlinkFormat")
    gwasfilename      <- paste0("GWASbetaStandard_allchrs", "sim", sim, ".txt")
    
    outnames_vec <- c()
    mafnames_vec <- c()
    for (ii in 1:3) {
      pop <- popvec[ii]
      base <- paste0(dirPlinkFormat, "/", pop, "chr", chr, "bedformat_numsnp", numsnp)
      outnames_vec <- c(outnames_vec, paste0(base, ".ld"))
    }
    
    Nvec <- TrainingNsam
    GWASbetafile <- file.path(dirSimuData, gwasfilename)
    
     summaryZoutput = Get_summary_3(GWASbetafile, Nvec, chr)
    summaryZ       <- summaryZoutput[[1]]
    Allkeepsnps    <- summaryZoutput[[2]]
    rm(summaryZoutput)
    
    plinkLD <- Gen_All_lddata_simpleTrans3(outnames_vec, Allkeepsnps, rcut = 0.03)
    
    # Generate Beta matrices
  
 savefile1 = paste0(output_sub_folder, penalty,"chr",chr,"_3pop","warmStart",warmStart,"sim",sim,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
    
      BetaMatrixS = Gen_One_BetaMatrix(savefile1, 3,1)[[1]]
      BetaMatrix3 = Gen_One_BetaMatrix(savefile1, 3,2)[[1]]
      BetaMatrixE = Gen_One_BetaMatrix(savefile1, 3,3)[[1]]

      mat = match(rownames(summaryZ),colnames(BetaMatrixS))
      BetaMatrixS = BetaMatrixS[,mat]

      mat = match(rownames(summaryZ),colnames(BetaMatrix3))
      BetaMatrix3 = BetaMatrix3[,mat]
      
      
      mat = match(rownames(summaryZ),colnames(BetaMatrixE))
      BetaMatrixE = BetaMatrixE[,mat]
      

      Start = 1
      for(i in 1:nrow(BetaMatrixS)){
        betavec = c(BetaMatrixS[i,],BetaMatrix3[i,],BetaMatrixE[i,])
        if(Start==1){
          AllBetaMatrix = betavec
          Start = 0
        }else{
          AllBetaMatrix = rbind(AllBetaMatrix, betavec)
        }
        rm(betavec)
      }
      P = nrow(summaryZ)
      
# Initialize an empty matrix to store the results
es_alphaMatrix <- matrix(0, nrow = 100, ncol = 8)

# Loop through each row of the matrices
for (i in 1:100) {
  
  # Extract the vectors for each population (row-wise)
  vec1 = BetaMatrixS[i, ]
  vec2 = BetaMatrix3[i, ]
  vec3 = BetaMatrixE[i, ]
  
  # Calculate the proportions for each of the 8 combinations
  pi_0 = mean(vec1 == 0 & vec2 == 0 & vec3 == 0)
  pi_1 = mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
  pi_2 = mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
  pi_3 = mean(vec1 != 0 & vec2 != 0 & vec3 == 0)
  pi_12 = mean(vec1 != 0 & vec2 != 0 & vec3 == 0)
  pi_13 = mean(vec1 != 0 & vec2 == 0 & vec3 != 0)
  pi_23 = mean(vec1 == 0 & vec2 != 0 & vec3 != 0)
  pi_123 = 1-(pi_0+pi_1+pi_2+pi_3+pi_12+pi_13+pi_23)
  
  # Store the results in the result matrix (row i)
  es_alphaMatrix[i, ] <- c(pi_0, pi_1, pi_2, pi_3, pi_12, pi_13, pi_23, pi_123)
}

    
    rm(BetaMatrix3, BetaMatrixS,BetaMatrixE)
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
  plinkLD = plinkLD,
  weight = 1.0
)
     
       
          SNPnames = rownames(summaryZ)
        P = length(SNPnames)
        inputs = output[[1]]
        #print(P)
        wBMatrix = apply(inputs,1,gen_sPostB2,P, num_alpha)
        rm(inputs)
        #rownames(wBMatrix1) = SNPnames
        Allrho_WBMatrix_list[[rr]] = wBMatrix
        rm(output, wBMatrix)
        gc()
      }
      rm(es_alphaMatrix,AllBetaMatrix)
      save(Allrho_WBMatrix_list,file=subtau_saveoutfile)
      rm(Allrho_WBMatrix_list)
    }
    rm(plinkLD, summaryZ, Allkeepsnps)
    gc()
  }
