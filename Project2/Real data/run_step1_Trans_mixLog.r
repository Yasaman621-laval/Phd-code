
library(dplyr)

  jset <-  as.numeric(args[1])
  inputs = read.table(paste0(file.rjobs,"input.txt"),as.is=T,header=T)

  runIndex = inputs[jset,2]
  chr = inputs[jset,3]
  print(chr)

  savename = "Trans2-realdata"
  penalty="RealmixLOG"
  filename = "run_step1_Trans_mixLog"
  NumL = 10
  subNumL = 10
  WeightN = weightN = 0


  
  Zscale = 1
  warmStart = 1
  singleStart = 1


 
  
  popvec = c("AFR","EAS","EUR")
  dirOutput = "/home/yatah3/projects/def-thchlava/yatah3/real-data/"
  output_sub_folder = paste0(dirOutput,savename,"/")
  system(paste0("mkdir -p ",output_sub_folder))

functionsfolder  <-"/home/yatah3/projects/def-thchlava/yatah3/real-data/Rfunction/"



  ordersequse_vec = c(1, 1)
  all_usedtrait = c("1,3","2,3")

  ordersequse = ordersequse_vec[runIndex]
  usedtrait = all_usedtrait[runIndex]

  usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
  usedtraitIndex = as.numeric(usedtraitsvec)


  mainindex = usedtraitIndex[1]
  

popuseY = popvec[mainindex]
  popuseT = popvec[3]

  library(SummaryLasso)
  savefile0 = NULL
  savef0 =  NULL

r2cut = 0.05

  NumIter = 1000

 source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))
    source(paste0(functionsfolder,"PRS_utility.r"))
    source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
    source(paste0(functionsfolder,"PlinkLD_transform.R"))
     source(paste0(functionsfolder,"Internal_Rfunctions.r"))
     
     

   #_______________________________________
#standardized Gwas colnames
#___________________________________________
#library(dplyr)
 
#out_dir = "/lustre06/project/6005709/yatah3/real-data/out/"
#GWASbetafile = paste0(out_dir,"final_clumped_crosspop.tsv")
    
 
#GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)

#colnames(GWASbeta) <- gsub("_AFR$", "_1", colnames(GWASbeta))
#colnames(GWASbeta) <- gsub("_EAS$", "_2", colnames(GWASbeta))
#colnames(GWASbeta) <- gsub("_EUR$", "_3", colnames(GWASbeta))


#write.table(GWASbeta,
     #  file = file.path(out_dir, "final_clumped_crosspop.tsv"),
        # sep = "\t",
        # quote = FALSE,
        # row.names = FALSE)

#_______________________________________________
   
   
out_dir = "/lustre06/project/6005709/yatah3/real-data/out/"
GWASbetafile = paste0(out_dir,"final_clumped_crosspop.tsv")
    
 
GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
   
   TrainingNsam = c(mean(GWASbeta$N_1),mean(GWASbeta$N_2),mean(GWASbeta$N_3))
   TrainingNsam <- round(TrainingNsam)
   Nvec = TrainingNsam[usedtraitIndex]
   
   summaryZoutput = Get_summaryZ(GWASbetafile, mainindex, Nvec, usedtraitsvec, chr)
    summaryZ = summaryZoutput[[1]]

  



  maxZ  =max(summaryZ$Z_3)



  low0_single = 0
    top0_single = maxZ

  
    for(ii in 1:2){
    Taufile = paste0(dirOutput,popvec[ii],"_",popvec[3],"Tau_info_v2.RData")

   if(!file.exists(Taufile)){
      AbsTauvec = matrix(,5)
       
GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
GWASbeta  = GWASbeta[,paste0("Z_",c(ii,3))]
        weights = apply(abs(GWASbeta),1,sum,na.rm=T)

  AbsTauvec = c(signif(quantile(weights,0)/10,3),
     signif(quantile(weights,0)/5,3),
      signif(quantile(weights,0)/2,3),
       signif(quantile(weights,0),3),
      signif(quantile(weights,0.1),3))
      save(AbsTauvec,file=Taufile)
    }
 
}

  
Taufile = paste0(dirOutput,popvec[mainindex],"_",popvec[3],"Tau_info_v2.RData")
load(file=Taufile)

  
   
 

r2cut=0.05


      if(Zscale==1){
        gwasfilename = GWASbetafile
      }
     outnames_vec = c()

      for(ii in 1:length(usedtraitIndex)){
        iiIndex = usedtraitIndex[ii]
        pop = popvec[iiIndex]
        outname = paste0(out_dir,"clumped_chr",chr,"_r2cut",r2cut,"_",pop)
        ld_outname = paste0(outname,".ld")
        outnames_vec = c(outnames_vec,ld_outname)
      }



   TransPEN_MixLog(Zscale = Zscale,usedtrait = usedtrait,TrainingNsam = TrainingNsam, dirSimuData = dirOutput, warmStart = warmStart,penalty = penalty,NumIter = NumIter,outnames_vec = outnames_vec,mainindex = mainindex,RupperVal = NULL, gwasfilename = gwasfilename, singleStart = singleStart, output_sub_folder_iter = output_sub_folder,AbsTauvec = AbsTauvec, savefile0 = savefile0, savef0 = savef0,
  chr = chr,
  weightN = WeightN,
  NumL = NumL,
  subNumL = subNumL,
  low0_single = 0,
  top0_single = top0_single , r2cut)



