rm(list = ls())

suppressPackageStartupMessages({
  library(SummaryLasso)
  library(data.table)
  library(dplyr)
  library(stringr)
})

# ============================================================
# 1) Job inputs (from SLURM array)
# ============================================================
args  <- commandArgs(trailingOnly = TRUE)
jset <- as.numeric(args[1])

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project1/step1/"
inputs <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE, as.is = TRUE)


sim           <- inputs[jset, 1]
chrIndex           <- inputs[jset, 2]
scenarioIndex <- inputs[jset, 3]

cat(">>> sim =", sim, " chr =", chrIndex, " scenario =", scenarioIndex, "\n")

# ============================================================
# 2) Global config
# ============================================================
nsim       <- 10
savename   <- "Trans2-new-3traits"
penalty    <- "RealmixLOG"
NumL       <- 10
subNumL    <- 10
WeightN    <- 0
Zscale     <- 1
rho_vec      <- c(seq(0,0.9,0.1),0.95)
K            <- 3
warmStart  <- 1
singleStart<- 1
NumIter    <- 1000

# population order (fixed)
popvec <- "EUR"
 
# ============================================================
# 3) Scenario selection + TrainingNsam parsing
#    (expects directories like .../SimuGenotype/sim_hsq..._train20000-20000-80000)
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/three_traits/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs)) stop("scenarioIndex exceeds available scenarios (", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)
cat(">>> Using scenario:", sc_name, "\n")

# parse TrainingNsam from folder name suffix after "_train"
train_str <- sub(".*_train", "", sc_name)            # "EUR50000"
train_num <- as.numeric(gsub("[^0-9]", "", train_str))  # remove non-digits ? "50000"

TrainingNsam    <- train_num

cat("    TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")
# ============================================================
# 4. General settings
# ============================================================
functions_folder <- "/lustre09/project/6005709/yatah3/simulation/Rfunction/"

dirSimuData  <- paste0(sc_dir, "/")
output_sub_folder <- file.path(sc_dir, "Trans2-new-3traits/")


output_sub_folder_uDensity <- file.path(sc_dir, savename, "/")
system(paste0("mkdir -p ", output_sub_folder_uDensity))

# ============================================================
# 5. Load variance of effect sizes from new summary .csv files
# ============================================================
output_sub_folder0 <- paste0("/lustre09/project/6005709/yatah3/simulation/SimuGenotype/three_traits/",sc_name,"/outputvisualonepop_sim",sim,"/")


var_beta <- numeric(3)

trait <- c("trait1", "trait2", "trait3")   

for (ii in 1:3) {

  finaloutputfit <- paste0(output_sub_folder0, trait[ii], "_sim", sim, ".csv")

  fit <- read.table(finaloutputfit, sep = ",", header = TRUE, as.is = TRUE)

  # corrected index
  var_beta[ii] <- fit$sig2_beta[1]
}


sigma2Kvec <- var_beta

cat("sigma2Kvec:", sigma2Kvec, "\n")


Taufile <- paste0(sc_dir,"/",popvec,"3traits_Tau_info_v2.RData")
if (!file.exists(Taufile)) stop("Taufile missing:", Taufile)
load(Taufile)  # loads AbsTauvec

# ============================================================
# 8. Load helper functions
# ============================================================
source(paste0(functions_folder, "AllPRS_Rfunctions.r"))
source(paste0(functions_folder, "PRS_utility.r"))
source(paste0(functions_folder, "Iterative_Rfunctions.r"))
source(paste0(functions_folder, "PlinkLD_transform.R"))

SharedPattern <- as.matrix(
  expand.grid(rep(list(c(0,1)), K))
)

  num_alpha = nrow(SharedPattern)





for (chr in chrIndex) {
  print(paste0(">>> chr=", chr))
  
  for (tt in 1:ncol(AbsTauvec)) {
    print(paste0(">>> Tau column=", tt))
    tauuse <- AbsTauvec[sim, tt]
    
    subtau_saveoutfile <- paste0(output_sub_folder_uDensity, penalty, "_chr", chr, 
                                 "_3traits", "_warmStart", warmStart, 
                                 "_sim", sim, "_Zscale", Zscale, "_singleStart", singleStart, 
                                 "_tauuse", tauuse, "_DensityU.RData")
    
    dirSimuDataSetChr <- file.path(dirBase, paste0("Chr", chr))
    dirPlinkFormat    <- file.path(dirSimuDataSetChr, "PlinkFormat")
    gwasfilename      <- paste0("GWASbetaStandard_allchrs_3traits_sim", sim, ".txt")
    
    outnames_vec <- c()
    mafnames_vec <- c()
  
      pop <- popvec
      base <- paste0(dirPlinkFormat, "/", pop, "chr",chr,"bedformat_numsnp400")
      outnames_vec <- c(outnames_vec, paste0(base, ".ld"))
    
    
    Nvec <- rep(TrainingNsam,3)
    GWASbetafile <- file.path(dirSimuData, gwasfilename)
    
    summaryZoutput <- Get_summary_3(GWASbetafile,Nvec,chr)
    summaryZ       <- summaryZoutput[[1]]
    Allkeepsnps    <- summaryZoutput[[2]]
    rm(summaryZoutput)
    
    plinkLD <- Gen_All_lddata_simpleTrans3traits(outnames_vec, Allkeepsnps, rcut = 0.03)
    
    # Generate Beta matrices
    savefile1 = paste0(output_sub_folder, penalty,"chr",chr,"_3traits","warmStart",warmStart,"sim",sim,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
    
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
