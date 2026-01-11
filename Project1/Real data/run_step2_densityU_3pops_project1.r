library(dplyr)

args <- commandArgs(trailingOnly = TRUE)

jset <-  as.numeric(args[1])
 
 file.rjobs = "/lustre06/project/6005709/yatah3/real-data/project1/"
 
 inputs = read.table(paste0(file.rjobs,"input.txt"),as.is=T,header=T)

 chrIndex <- inputs[jset, 1]

  
  savename = "new_densityU_3pops_project1"

 penalty="RealmixLOG"


  warmStart = 1
  K  = 3
  
  savename0 = "mixer"
  
 
  dirOutput = "/home/yatah3/projects/def-thchlava/yatah3/real-data/project1/"

  
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

  sigma2Kvec = var_beta


  Zscale = 1
  library(SummaryLasso)
library(gtools)

  dirtemp = Sys.getenv('SLURM_TMPDIR')
  folder_predi = paste0(dirtemp,"/")

out_dir = "/lustre06/project/6005709/yatah3/real-data/out/"
GWASbetafile = paste0(out_dir,"final_clumped_crosspop.tsv")
    
 
GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)


   TrainingNsam = c(mean(GWASbeta$N_1),mean(GWASbeta$N_2),mean(GWASbeta$N_3))
   TrainingNsam <- round(TrainingNsam)
  
  
functionsfolder  <- "/lustre06/project/6005709/yatah3/real-data/Rfunction/"

  source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))
  source(paste0(functionsfolder,"PRS_utility.r"))
  source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
  source(paste0(functionsfolder,"PlinkLD_transform.R"))





 
 
 
  output_sub_folder = "/lustre06/project/6005709/yatah3/real-data/project1/Trans2-new-3-out/"

 
output_sub_folder_uDensity = paste0("/lustre06/project/6005709/yatah3/real-data/project1/", savename, "/")
system(paste0("mkdir -p '", output_sub_folder_uDensity, "'"))


 

 



    all_usedtrait = "1,2,3"
  usedtrait = all_usedtrait
  usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
  usedtraitIndex = as.numeric(usedtraitsvec)

 




  rho_vec = c(seq(0,0.9,0.1),0.95)
  groupnum = 3
  singleStart = 1

  SharedPattern = permutations(n=2,r=K,v=c(0,1),repeats.allowed=T)
  num_alpha = nrow(SharedPattern)


  Taufile =paste0(dirOutput,popvec[1],"_",popvec[2], "-", popvec[3],"Tau_info_v2.RData")
  load(file=Taufile)




chr = chrIndex
  print(paste0("chr=", chr))
    for(tt in 1:length(AbsTauvec)){
     print(paste0("tt=", tt))
      tauuse = AbsTauvec[tt]
      subtau_saveoutfile =  paste0(output_sub_folder_uDensity, penalty, "chr", chr, 
                       "usedtrait_1,2,3", "warmStart", warmStart, 
                        "Zscale", Zscale, "singleStart", singleStart, 
                       "tauuse", tauuse, "_DensityU.RData")
                       
  
   rcut=0.05

      
   outnames_vec = c()

      for(ii in 1:length(usedtraitIndex)){
        iiIndex = usedtraitIndex[ii]
        pop = popvec[iiIndex]
        outname = paste0(out_dir,"clumped_chr",chr,"_r2cut",rcut,"_",pop)
        ld_outname = paste0(outname,".ld")
        outnames_vec = c(outnames_vec,ld_outname)
      }


    Nvec = TrainingNsam

    summaryZoutput = Get_summary_3(GWASbetafile,Nvec, usedtraitsvec, chr)
    summaryZ = summaryZoutput[[1]]
    Allkeepsnps = summaryZoutput[[2]]
    rm(summaryZoutput)

    plinkLD = Gen_All_lddata_simpleTrans3(outnames_vec, Allkeepsnps, rcut = 0.05)
   
  


      savefile1 = paste0(output_sub_folder,penalty,"chr",chr,"_3traits","warmStart",warmStart,"Zscale",Zscale,"tauuse",tauuse,"singleStart",singleStart,".RData")
 

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
      cat("Saving to", subtau_saveoutfile, "\n")
      rm(Allrho_WBMatrix_list)
    }
