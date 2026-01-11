#!/usr/bin/env Rscript
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
chr           <- inputs[jset, 2]
scenarioIndex <- inputs[jset, 3]

cat(">>> sim =", sim, " chr =", chr, " scenario =", scenarioIndex, "\n")

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
# 4) Paths
# ============================================================
dirSimuDataSet <- paste0(sc_dir, "/")              # scenario root
output_sub_folder <- file.path(sc_dir, savename, "/")   # results folder per scenario
system(paste0("mkdir -p ", shQuote(output_sub_folder)))

functions_folder <- "/lustre09/project/6005709/yatah3/simulation/Rfunction/"
 
Taufile = paste0(sc_dir,"/",popvec,"3traits_Tau_info_v2.RData")

   if(!file.exists(Taufile)){
      AbsTauvec = matrix(,nsim,5)
      for(sim in 1:nsim){
        GWASbetafile = file.path(sc_dir, paste0("GWASbetaStandard_allchrs_3traits_sim", sim, ".txt"))
        GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
        GWASbeta  = GWASbeta[,paste0("b",c(1,2,3))]
        weights = apply(abs(GWASbeta),1,sum,na.rm=T)

        AbsTauvec[sim,1] = signif(quantile(weights,0)/10,3)
        AbsTauvec[sim,2] = signif(quantile(weights,0)/5,3)
        AbsTauvec[sim,3] = signif(quantile(weights,0)/2,3)
        AbsTauvec[sim,4] = signif(quantile(weights,0),3)
        AbsTauvec[sim,5] = signif(quantile(weights,0.1),3)
        rm(GWASbeta)
      }

  save(AbsTauvec, file = Taufile)
}


load(Taufile)  # provides AbsTauvec

# ============================================================
# 7) Load helper functions & compute summary Z for chosen chr
# ============================================================
source(file.path(functions_folder, "AllPRS_Rfunctions.r"))
source(file.path(functions_folder, "PRS_utility.r"))
source(file.path(functions_folder, "Iterative_Rfunctions.r"))
source(file.path(functions_folder, "PlinkLD_transform.R"))
source(file.path(functions_folder, "Internal_Rfunctions.r")) 

GWASbetafile = file.path(sc_dir, paste0("GWASbetaStandard_allchrs_3traits_sim", sim, ".txt"))
   summaryZoutput = Get_summary_3(GWASbetafile, TrainingNsam, chr)
    summaryZ = summaryZoutput[[1]]



  low0_single = 0

   maxZ1 = max(summaryZ$Zobs1)
  maxZ2  =max(summaryZ$Zobs2)
 maxZ3  =max(summaryZ$Zobs3)

  RupperVal1 = max(c((maxZ1/sqrt(TrainingNsam))*2))
  RupperVal2 = max(c((maxZ2/sqrt(TrainingNsam))*2))
RupperVal3 =  max(c((maxZ3/sqrt(TrainingNsam))*2))

  low0_single = 0
  top0_single1 = maxZ1
   top0_single2 = maxZ2
  top0_single3 = maxZ3
  
  
  
  
  
   
    dirSimuDataSetChr = paste0(dirBase,"Chr",chr,"/")
    dirPlinkFormat = paste0(dirSimuDataSetChr,"PlinkFormat/")


      outnames_vec = c()
      mafnames_vec = c()
        pop = popvec
 
        outname = paste0(dirPlinkFormat,pop,"chr",chr,"bedformat_numsnp400")
        ld_outname = paste0(outname,".ld")
        outnames_vec = c(outnames_vec,ld_outname)
        maffile =  paste0(dirPlinkFormat,pop,"chr",chr,"bedformat.frq")
        mafnames_vec = c(mafnames_vec,maffile)
   


      TransPEN_MixLog3traits( Zscale        = Zscale,
  usedtrait     = usedtrait,           
  TrainingNsam  = TrainingNsam,
  dirSimuData   = dirSimuDataSet,        
  warmStart     = warmStart,
  penalty       = penalty,
  NumIter       = NumIter,
  outnames_vec  = outnames_vec,        
  mainindex     = mainindex,         
  RupperVal     = RupperVal,
  gwasfilename  = GWASbetafile,
  mafnames_vec  = mafnames_vec,        
  singleStart   = singleStart,
  output_sub_folder = output_sub_folder,
  AbsTauvec     = AbsTauvec,
  sim           = sim,
  savefile0     = savefile0,
  savef0        = savef0,
  chr           = chr,
  WeightN       = WeightN,
  NumL          = NumL,
  subNumL       = subNumL,
  low0_single   = low0_single,
   top0_single1 = top0_single1,
   top0_single2 = top0_single2,
   top0_single3 = top0_single3 )
  
  


