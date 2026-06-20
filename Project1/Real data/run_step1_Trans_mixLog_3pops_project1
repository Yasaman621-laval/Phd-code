



args <- commandArgs(trailingOnly = TRUE)

  jset <-  as.numeric(args[1])
  
  file.rjobs = "/lustre06/project/6005709/yatah3/real-data/project1/"
  
  inputs = read.table(paste0(file.rjobs,"input.txt"),as.is=T,header=T)

  chr = inputs[jset,1]

  print(chr)

  savename = "Trans2-new-3-out"
  penalty="RealmixLOG"
  filename = "run_step1_Trans_mixLog"
  NumL = 10
  subNumL = 10
  WeightN = weightN = 0

  
  Zscale = 1
  warmStart = 1
  singleStart = 1


 
  
  popvec = c("AFR","EAS","EUR")
  dirOutput = "/lustre06/project/6005709/yatah3/real-data/project1/"
  output_sub_folder = paste0(dirOutput,savename,"/")
  system(paste0("mkdir -p ",output_sub_folder))



  usedtrait = "1,2,3"

  usedtraitsvec = unlist(strsplit(as.character(usedtrait),split=","))
  usedtraitIndex = as.numeric(usedtraitsvec)



  library(SummaryLasso)
  savefile0 = NULL
  savef0 =  NULL



  NumIter = 1000
  
 out_dir = "/lustre06/project/6005709/yatah3/real-data/out/"
GWASbetafile = paste0(out_dir,"final_clumped_crosspop.tsv")
    
 
GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
   
   TrainingNsam = c(mean(GWASbeta$N_1),mean(GWASbeta$N_2),mean(GWASbeta$N_3))
   TrainingNsam <- round(TrainingNsam)
 
  N3 = TrainingNsam[3]


functionsfolder  <- "/lustre06/project/6005709/yatah3/real-data/Rfunction/"




    Taufile = paste0(dirOutput,popvec[1],"_",popvec[2], "-", popvec[3],"Tau_info_v2.RData")

    if(!file.exists(Taufile)){
      AbsTauvec = matrix(,5)
        GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
        GWASbeta  = GWASbeta[,paste0("beta_",c(1,2,3))]
        weights = apply(abs(GWASbeta),1,sum,na.rm=T)

    AbsTauvec = c(signif(quantile(weights,0)/10,3),
     signif(quantile(weights,0)/5,3),
      signif(quantile(weights,0)/2,3),
       signif(quantile(weights,0),3),
      signif(quantile(weights,0.1),3))
      save(AbsTauvec,file=Taufile)
    }
 


load(file=Taufile)

 source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))
    source(paste0(functionsfolder,"PRS_utility.r"))
    source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
    source(paste0(functionsfolder,"PlinkLD_transform.R"))
     source(paste0(functionsfolder,"Internal_Rfunctions.r"))
     

   summaryZoutput = Get_summary_3(GWASbetafile, TrainingNsam, usedtraitsvec, chr)
    summaryZ = summaryZoutput[[1]]


  maxZ1 = max(summaryZ$Z_1)
  maxZ2 = max(summaryZ$Z_2)
  maxZ3 = max(summaryZ$Z_3)

 

  low0_single = 0
  
  top0_single1 = maxZ1
   top0_single2 = maxZ2
  top0_single3 = maxZ3
  
  
   rcut=0.05

      
   outnames_vec = c()

      for(ii in 1:length(usedtraitIndex)){
        iiIndex = usedtraitIndex[ii]
        pop = popvec[iiIndex]
        outname = paste0(out_dir,"clumped_chr",chr,"_r2cut",rcut,"_",pop)
        ld_outname = paste0(outname,".ld")
        outnames_vec = c(outnames_vec,ld_outname)
      }



      TransPEN_MixLog3(Zscale, TrainingNsam, dirSimuData = dirOutput, warmStart, penalty, NumIter, outnames_vec, mainindex, GWASbetafile, singleStart, output_sub_folder, AbsTauvec, savefile0, savef0, chr, WeightN, NumL = NumL, subNumL = subNumL, low0_single = 0, top0_single1,top0_single2,top0_single3,rcut)
  
  


