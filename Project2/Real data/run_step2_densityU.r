library(dplyr)
jset <-  as.numeric(args[1])
  load(paste0(file.rjobs,"input.RData"))
  inputs = Alltemp

  runIndex =inputs[jset,1]
  chrIndex = inputs[jset,2]
  

  savename = "densityU_realdata"

 
  testset = 2


  warmStart = 1
  K  = 2
  
  
  dirOutput = "/home/yatah3/projects/def-thchlava/yatah3/real-data/"

  
 output_sub_folder0 = "/home/yatah3/projects/def-thchlava/yatah3/real-data/runR-mixer/"
  
  popvec = c("AFR","EAS","EUR")
var_beta <- NULL  # Initialize var_beta if not already initialized




for (ii in 1:3) {
  finaloutputfit <- paste0(output_sub_folder0,"sig2_beta_per_chr_", popvec[ii],".csv")
fit <- read.table(
  file = finaloutputfit,
  sep = ",",          # ? use comma, not tab
  header = TRUE
)
var = fit[nrow(fit),2]
var_beta <- rbind(var_beta, var)
}

var_beta = as.vector(unlist(var_beta))

  Zscale = 1
  library(SummaryLasso)

  dirtemp = Sys.getenv('SLURM_TMPDIR')
  folder_predi = paste0(dirtemp,"/")
  ordersequse_vec = c(1, 1)
  all_usedtrait = c("1,3","2,3")

out_dir = "/lustre06/project/6005709/yatah3/real-data/out/"
GWASbetafile = paste0(out_dir,"final_clumped_crosspop.tsv")
    
 
GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)


   TrainingNsam = c(mean(GWASbeta$N_1),mean(GWASbeta$N_2),mean(GWASbeta$N_3))
   TrainingNsam <- round(TrainingNsam)
  
  
  
  source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))
  source(paste0(functionsfolder,"PRS_utility.r"))
  source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
  source(paste0(functionsfolder,"PlinkLD_transform.R"))



  ordersequse_vec = c(1, 1)
  all_usedtrait = c("1,3","2,3")

  ordersequse = ordersequse_vec[runIndex]
  usedtrait = all_usedtrait[runIndex]

  usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
  usedtraitIndex = as.numeric(usedtraitsvec)


  mainindex = usedtraitIndex[1]

 
 
 
  output_sub_folder = "/home/yatah3/projects/def-thchlava/yatah3/real-data/Trans2-realdata/"
 
 
 

  output_sub_folder_uDensity = paste0(dirOutput,savename,"/")
  system(paste0("mkdir -p ",output_sub_folder_uDensity))




  if(runIndex==1){
  sigma2Kvec = var_beta[c(1,3)]
  }else{
  sigma2Kvec = var_beta[c(2,3)]
  }

#sigma2Kvec = var_beta[c(1,2)]





  N = TrainingNsam[3]

  
  rho_vec = c(seq(0,0.9,0.1),0.95)
  groupnum = 2
  singleStart = 1

  SharedPattern = permutations(n=2,r=K,v=c(0,1),repeats.allowed=T)
  num_alpha = nrow(SharedPattern)



Taufile = paste0(dirOutput,popvec[mainindex],"_",popvec[3],"Tau_info_v2.RData")
load(file=Taufile)


r2cut = 0.05

for (chr in chrIndex_list[[chrIndex]]) {
  print(paste0("chr=", chr))
    for(tt in 1:length(AbsTauvec)){
     print(paste0("tt=", tt))
      tauuse = AbsTauvec[tt]
      subtau_saveoutfile =  paste0(output_sub_folder_uDensity, penalty, "chr", chr, 
                       "usedtrait", usedtrait, "warmStart", warmStart, 
                        "Zscale", Zscale, "singleStart", singleStart, 
                       "tauuse", tauuse, "_DensityU.RData")
                       
     #if (file.exists(subtau_saveoutfile) && file.info(subtau_saveoutfile)$size > 100 * 1024 * 1024) {next}

       dirSimuDataSetChr = paste0(dirOutput,"Chr",chr,"/")
    dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")
    
    gwasfilename = paste0(out_dir,"final_clumped_crosspop.tsv")

    outnames_vec = c()
   
    for(ii in 1:length(usedtraitIndex)){
        iiIndex = usedtraitIndex[ii]
        pop = popvec[iiIndex]
        outname = paste0(out_dir,"clumped_chr",chr,"_r2cut",r2cut,"_",pop)
        ld_outname = paste0(outname,".ld")
        outnames_vec = c(outnames_vec,ld_outname)
      }

    Nvec = TrainingNsam[usedtraitIndex]


    summaryZoutput = Get_summaryZ(gwasfilename, mainindex, Nvec, usedtraitsvec, chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
    rm(summaryZoutput)

    plinkLD = Gen_All_lddata_simpleTrans(outnames_vec, Allkeepsnps, rcut = 0.05)
   
   

      savefile1 = paste0(output_sub_folder, penalty,"chr",chr,"usedtrait",usedtrait,"warmStart",warmStart,"Zscale",Zscale,"tauuse",tauuse,".RData")



      BetaMatrixS = Gen_One_BetaMatrix(savefile1, 2,1)[[1]]
      BetaMatrix3 = Gen_One_BetaMatrix(savefile1, 2,2)[[1]]

      mat = match(rownames(summaryZ),colnames(BetaMatrixS))
      BetaMatrixS = BetaMatrixS[,mat]

      mat = match(rownames(summaryZ),colnames(BetaMatrix3))
      BetaMatrix3 = BetaMatrix3[,mat]

      Start = 1
      for(i in 1:nrow(BetaMatrixS)){
        betavec = c(BetaMatrixS[i,],BetaMatrix3[i,])
        if(Start==1){
          AllBetaMatrix = betavec
          Start = 0
        }else{
          AllBetaMatrix = rbind(AllBetaMatrix, betavec)
        }
        rm(betavec)
      }
      P = nrow(summaryZ)
      
# Initialize the results matrix
es_alphaMatrix <- matrix(nrow =nrow(AllBetaMatrix), ncol = 4)

# Loop through each row of AllBeta
for (i in 1:nrow(AllBetaMatrix)) {
  # Extract non-zero elements from AllBeta[i,]
  non_zero_elements <- AllBetaMatrix[i,][AllBetaMatrix[i,] != 0]

  if (length(non_zero_elements) == 0) {
    # If there are no non-zero elements, return a vector of zeros
    result2 <- c(1, 0, 0, 0)
  } else {
    # Number of non-zero elements required
    num_non_zero_required <- length(non_zero_elements)
    
    # Calculate half the number of non-zero elements for each column
    half_non_zero <- floor(num_non_zero_required /length(Nvec))
    
    # Initialize reshaped_beta with zeros
    reshaped_beta <- matrix(0, nrow = ncol(BetaMatrixS), ncol = length(Nvec))
    
    # Set the non-zero elements in the reshaped_beta
    reshaped_beta[1:half_non_zero, 1] <- non_zero_elements[1:half_non_zero]
    reshaped_beta[1:(num_non_zero_required - half_non_zero), length(Nvec)] <- non_zero_elements[(half_non_zero + 1):num_non_zero_required]
    
    # Ensure reshaped_beta has the same number of non-zero values
    vec1 <- reshaped_beta[, 1]
    vec2 <- reshaped_beta[, 2]
    
    # Compute the required means
    mean_vec1_eq_0_vec2_eq_0 <- mean(vec1 == 0 & vec2 == 0)
    mean_vec1_eq_0_vec2_neq_0 <- mean(vec1 == 0 & vec2 != 0)
    mean_vec1_neq_0_vec2_eq_0 <- mean(vec1 != 0 & vec2 == 0)
    mean_vec1_neq_0_vec2_neq_0 <- mean(vec1 != 0 & vec2 != 0)
    
    # Combine the results into a vector
    result2 <- c(mean_vec1_eq_0_vec2_eq_0, mean_vec1_eq_0_vec2_neq_0, mean_vec1_neq_0_vec2_eq_0, mean_vec1_neq_0_vec2_neq_0)
  }
  
  # Assign the result to the results matrix
  es_alphaMatrix[i, ] <- result2
}
rm(BetaMatrix3,BetaMatrixS)
gc()
 
      Allrho_WBMatrix_list = list()
      for(rr in 1:length(rho_vec)){
        rho = rho_vec[rr]
        rhoMat = matrix(rho,K,K)
        sigma2K_allAlpha_List = Create_sigma2K_allAlpha_List(K,sigma2Kvec,rhoMat)

        output = transConditionalU(summaryZ = summaryZ, Nvec = Nvec, JointBmatrix = AllBetaMatrix, Zcov = diag(K), sigma2K_allAlpha_List=sigma2K_allAlpha_List, es_alphaMatrix= es_alphaMatrix, plinkLD=plinkLD)
        
        
     SNPnames = rownames(summaryZ)
        P = length(SNPnames)
        inputs = cbind(output[[1]],output[[4]])
        wBMatrix1 = apply(inputs,1,gen_sPostB,P, num_alpha)
        rm(inputs)
        rownames(wBMatrix1) = SNPnames
        inputs = cbind(output[[1]],output[[5]])
        wBMatrix2 = apply(inputs,1,gen_sPostB,P, num_alpha)
        rownames(wBMatrix2) = SNPnames
        rm(inputs)
        wBoutput = list()
        wBoutput[[1]] = wBMatrix1
        wBoutput[[2]] = wBMatrix2
        rm(output, wBMatrix1,wBMatrix2)
        gc()
        Allrho_WBMatrix_list[[rr]] = wBoutput
        rm(wBoutput)
      }
      rm(es_alphaMatrix,AllBetaMatrix)
      save(Allrho_WBMatrix_list,file=subtau_saveoutfile)
      rm(Allrho_WBMatrix_list)
    }
    rm(plinkLD, summaryZ, Allkeepsnps)
    gc()
  }
