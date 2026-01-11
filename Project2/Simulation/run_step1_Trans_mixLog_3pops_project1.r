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

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project2/step1/"
inputs <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE, as.is = TRUE)


sim           <- inputs[jset, 1]
runIndex      <- inputs[jset, 2]  # 1 => AFR+EUR, 2 => EAS+EUR
chr           <- inputs[jset, 3]
scenarioIndex <- inputs[jset, 4]

cat(">>> sim =", sim, " runIndex =", runIndex, " chr =", chr, " scenario =", scenarioIndex, "\n")

# ============================================================
# 2) Global config
# ============================================================
nsim       <- 10
savename   <- "Trans2-new-3"
penalty    <- "RealmixLOG"
NumL       <- 10
subNumL    <- 10
WeightN    <- 0
Zscale     <- 1
warmStart  <- 1
singleStart<- 1
NumIter    <- 2000

# population order (fixed)
popvec <- c("AFR","EAS","EUR")   # 1,2 are non-EUR; 3 is EUR

# ============================================================
# 3) Scenario selection + TrainingNsam parsing
#    (expects directories like .../SimuGenotype/sim_hsq..._train20000-20000-80000)
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs)) stop("scenarioIndex exceeds available scenarios (", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)
cat(">>> Using scenario:", sc_name, "\n")

# parse TrainingNsam from folder name suffix after "_train"
train_str       <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)
nums            <- as.numeric(sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x))))
TrainingNsam    <- nums

cat("    TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")

# ============================================================
# 4) Paths
# ============================================================
dirSimuDataSet <- paste0(sc_dir, "/")              # scenario root
output_sub_folder <- file.path(sc_dir, savename, "/")   # results folder per scenario
system(paste0("mkdir -p ", shQuote(output_sub_folder)))

functions_folder <- "/lustre09/project/6005709/yatah3/simulation/Rfunction/"

# ============================================================
# 5) Pick pair indices: (AFR, EUR) or (EAS, EUR) by runIndex
#    runIndex = 1 -> "1,3" (AFR,EUR)
#    runIndex = 2 -> "2,3" (EAS,EUR)
# ============================================================
all_usedtrait <- c("1,3", "2,3")
if (!(runIndex %in% 1:2)) stop("runIndex must be 1 or 2 for pairwise EUR runs.")

usedtrait     <- all_usedtrait[runIndex]
usedtraitsvec <- as.numeric(unlist(strsplit(usedtrait, ",")))  # length 2
mainindex     <- usedtraitsvec[1]  # AFR (1) or EAS (2) is the "main" non-EUR pop
eurindex      <- usedtraitsvec[2]  # should be 3

pop_main <- popvec[mainindex]  # "AFR" or "EAS"
pop_eur  <- popvec[eurindex]   # "EUR"

cat(">>> Pair =", paste(pop_main, pop_eur, sep = "+"),
    " (usedtrait:", usedtrait, ") mainindex =", mainindex, "\n")

# ============================================================
# 6) Tau grid for THIS PAIR (main, EUR)
#    File matches original style: "<main>_<EUR>Tau_info_v2.RData"
# ============================================================
 
Taufile = paste0(sc_dir,"/",popvec[mainindex],"_",popvec[eurindex],"Tau_info_v2.RData")

   if(!file.exists(Taufile)){
      AbsTauvec = matrix(,nsim,5)
      for(s in 1:nsim){
        GWASbetafile = file.path(sc_dir, paste0("GWASbetaStandard_allchrssim", s, ".txt"))
        GWASbeta = read.table(file=GWASbetafile, sep = "\t", header=T, as.is=T)
        GWASbeta  = GWASbeta[,paste0("b",c(mainindex,eurindex))]
        weights = apply(abs(GWASbeta),1,sum,na.rm=T)

        AbsTauvec[s,1] = signif(quantile(weights,0)/10,3)
        AbsTauvec[s,2] = signif(quantile(weights,0)/5,3)
        AbsTauvec[s,3] = signif(quantile(weights,0)/2,3)
        AbsTauvec[s,4] = signif(quantile(weights,0),3)
        AbsTauvec[s,5] = signif(quantile(weights,0.1),3)
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

GWASbetafile <- file.path(sc_dir, paste0("GWASbetaStandard_allchrssim", sim, ".txt"))
if (!file.exists(GWASbetafile)) stop("Missing GWAS file: ", GWASbetafile)


 Nvec = TrainingNsam[usedtraitsvec]

summaryZoutput <- Get_summaryZ(GWASbetafile,mainindex, Nvec, usedtraitsvec, chr)
summaryZ <- summaryZoutput[[1]]

# Map mainindex to the right Zobs column (1->Zobs1, 2->Zobs2, 3->Zobs3)
Zobs_col <- paste0("Zobs", mainindex)

maxZ_main <- max(summaryZ[[Zobs_col]], na.rm = TRUE)

N_main <- TrainingNsam[mainindex]
RupperVal <- max((maxZ_main / sqrt(N_main)) * 2)

low0_single <- 0
top0_single <- maxZ_main


# ============================================================
# 8) Prepare PLINK LD / MAF inputs ONLY for the chosen pair
# ============================================================
dirSimuDataSetChr <- file.path(sc_dir,"/", paste0("Chr", chr))
dirPlinkFormat    <- file.path(dirSimuDataSetChr, "PlinkFormat")

numsnp <- 400
pair_indices <- usedtraitsvec  # e.g., c(1,3) or c(2,3)

outnames_vec <- character(0)
mafnames_vec <- character(0)
for (ii in pair_indices) {
  pop <- popvec[ii]
  base <- file.path(dirPlinkFormat, paste0(pop, "chr", chr, "bedformat_numsnp", numsnp))
  outnames_vec <- c(outnames_vec, paste0(base, ".ld"))
  mafnames_vec <- c(mafnames_vec, file.path(dirPlinkFormat, paste0(pop, "chr", chr, "bedformat.frq")))
}


# ============================================================
# 9) Run TWO-pop Trans-ethnic MixLog
#    NOTE: this calls the 2-pop function, not the 3-pop version.
# ============================================================
savefile0 <- NULL
savef0    <- NULL

cat(">>> Running TransPEN_MixLog for pair ", pop_main, "+", pop_eur, " on chr ", chr, " ...\n", sep = "")
TransPEN_MixLog(
  Zscale        = Zscale,
  usedtrait     = usedtrait,             # "1,3" or "2,3"
  TrainingNsam  = TrainingNsam,
  dirSimuData   = dirSimuDataSet,        # scenario root
  warmStart     = warmStart,
  penalty       = penalty,
  NumIter       = NumIter,
  outnames_vec  = outnames_vec,          # only 2 LD files
  mainindex     = mainindex,             # 1 for AFR or 2 for EAS
  RupperVal     = RupperVal,
  gwasfilename  = GWASbetafile,
  mafnames_vec  = mafnames_vec,          # only 2 MAF files
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
  top0_single   = top0_single
)

cat(">>> DONE: sim", sim, "runIndex", runIndex, "chr", chr, "scenario", scenarioIndex, "\n")
