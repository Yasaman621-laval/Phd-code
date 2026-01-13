rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing SLURM_ARRAY_TASK_ID")

jset <- as.numeric(args[1])

input_file <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/input_scenario.txt"

inputs <- read.table(input_file, header = TRUE, as.is = TRUE)

scenarioIndex <- inputs[jset, 1]

cat(" scenario =", scenarioIndex, "\n")

savename <- "CSx-meta"

dirOutput <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/"

path_to_python <- "/lustre09/project/6005709/yatah3/simulation/project2/PRs-csx/PRScsx/"

popvec <- c("AFR", "EAS", "EUR")

# ============================================================
# 3. Load summary stats for sample sizes
# ============================================================
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) stop("No scenario directories found under: ", dirBase)
if (scenarioIndex > length(scenario_dirs)) stop("scenarioIndex exceeds available scenarios (", length(scenario_dirs), ")")

sc_dir  <- scenario_dirs[scenarioIndex]
dirSimuDataSet = paste0(sc_dir,"/")
output_sub_folder <- paste0(sc_dir,"/", savename, "/")

sc_name <- basename(sc_dir)
cat(">>> Using scenario:", sc_name, "\n")

# parse TrainingNsam from folder name suffix after "_train"
train_str       <- sub(".*_train", "", sc_name)
train_str_clean <- gsub("e\\+0", "e", train_str)
nums            <- as.numeric(sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x))))
TrainingNsam    <- nums

cat("    TrainingNsam =", paste(TrainingNsam, collapse = ", "), "\n")
# ============================================================
# 4. phi values (ttIndex maps here)
# ============================================================

phivec <- c(1e-6, 1e-4, 1e-2, 1)

PV = ""
phenoIndex = 3

nsim = 10

functionsfolder = "/lustre09/project/6005709/yatah3/Rfunction_PRS_CSX/"
source(paste0(functionsfolder,"PRS_utility.r"))
source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))

dirtemp = Sys.getenv('SLURM_TMPDIR')
folder_predi = paste0(dirtemp,"/")

R2output = paste0(sc_dir,"/",savename,"/R2output/")
system(paste0("mkdir -p ",R2output))
Scoreoutput = paste0(folder_predi,savename,"/Scoreoutput/")
system(paste0("mkdir -p ",Scoreoutput))

popvec = c("AFR","EAS","EUR")

refAllelefileInput = paste0(dirBase,"refereceAlleleSNPs.txt")
if(savename == "CSx"){
  refAllelefileInput = NULL
}

TrainIndex = 3
popTrain = popvec[TrainIndex]
# ============================================================
# ?? BEGIN MAIN SCRIPT (unchanged except for MAF mapping lines)
# ============================================================
############################################################
##  MASTER POST-PRSCSx PARSING PIPELINE (corrected)
##  - Handles EUR+AFR and EUR+EAS
##  - Handles NON-META and META outputs
##  - Uses unified reference SNP mapping (rawMAF)
############################################################

for (sim in 1:nsim) {

  ## Load unified SNP mapping (ALWAYS used for matching)
  rawMAF = read.table(
    file = paste0(dirBase, "matchSNP_allPops.txt"),
    as.is = TRUE, header = TRUE
  )

  for (iiIndex in 1:3) {

    pop = popuseY = popvec[iiIndex]   # AFR or EAS (EUR is skipped)
    if (pop == "EUR") next            # EUR is always the TRAINING population

    ## Select evaluation MAF columns (used later for R2)
    MAFcol = paste0("MAF_", pop)
    NCHcol = paste0("NCHROBS_", pop)

    MAFdata = data.frame(
      CHR     = rawMAF$CHR,
      SNP     = rawMAF$SNP,
      A1      = rawMAF$A1,
      A2      = rawMAF$A2,
      OA      = rawMAF$OA,
      MAF     = rawMAF[[MAFcol]],
      NCHROBS = rawMAF[[NCHcol]],
      stringsAsFactors = FALSE
    )

    for (phiuse in phivec) {

      base_prefix <- paste0("PRSCSx_", popTrain, "_", pop,
                            "_phi", phiuse, "_sim", sim)

      # Fix: require simX_ so sim1 does NOT match sim10
      pattern_sim <- paste0(base_prefix, "_.*_chr[0-9]+\\.txt$")

      finishedfiles_all <- list.files(
        path = output_sub_folder,
        pattern = pattern_sim,
        full.names = FALSE
      )

      if (length(finishedfiles_all) == 0) next

      files_nonmeta <- finishedfiles_all[!grepl("_META_", finishedfiles_all)]

      ############################################################
      ##  A) NON-META BRANCH  (EUR file + AFR/EAS file per chr)
      ############################################################
      if (length(files_nonmeta) > 0) {

        finishedchr <- c()
        for (fn in files_nonmeta) {
          parts   <- strsplit(fn, "_")[[1]]
          chrpart <- parts[grep("chr", parts)]
          chrnum  <- as.numeric(gsub("chr|\\.txt","", chrpart))
          finishedchr <- c(finishedchr, chrnum)
        }

        test <- table(finishedchr)

        ## Expect 2 files per chromosome (EUR + target)
        if (length(which(test == 2)) == 22) {

          savecore <- paste0(popTrain, "_", pop,
                             "sim", sim, "phiuse", phiuse)

         
  savef2_list = list()
            savef2_list[[1]] =  paste0(R2output, savecore,"Train_results.RData")
            savef2_list[[2]] =  paste0(R2output, savecore,"Test_results.RData")
      
     AllFinish = 1
            for(ll in 1:length(savef2_list)){
              savef2 = savef2_list[[ll]]
              if(file.exists(savef2)){
                load(savef2)
                if(Finish==0){
                  AllFinish = 0
                  break
                }
              }else{
                AllFinish = 0
                break
              }
            }
            
     if(AllFinish==0){
            output_sub_folder_predi <- paste0(folder_predi, savecore, "/")
            system(paste0("mkdir -p ", output_sub_folder_predi))

            for (chr in 1:22) {

              chrname = paste0("chr", chr, ".txt")
              idx     = grep(chrname, files_nonmeta)
              usefiles = files_nonmeta[idx]

              if (length(usefiles) != 2)
                stop("Expected 2 NON-META files for chr=", chr,
                     " but found ", length(usefiles))

              extract_pop <- function(filename) {
                # captures pop between _simX_ and _pst_eff
                sub(".*_sim[0-9]+_([^_]+)_pst_eff.*", "\\1", filename)
              }

              ## Identify which file is EUR and which is target AFR/EAS
              postpops <- sapply(usefiles, extract_pop)

              grepTrain <- which(postpops == popTrain)

              if (length(grepTrain) != 1)
                stop("Could not find exactly 1 EUR posterior file for chr = ", chr)

              grepTest <- setdiff(1:2, grepTrain)

              if (length(grepTest) != 1)
                stop("Could not find exactly 1 target-pop posterior file for chr = ", chr)

              fileTrain = paste0(output_sub_folder, usefiles[grepTrain])
              fileTest  = paste0(output_sub_folder, usefiles[grepTest])

              betaTrain <- read.table(fileTrain, as.is=TRUE)
              betaTest  <- read.table(fileTest,  as.is=TRUE)

              ## ================
              ## Correct mapping
              ## ================
              matTrain <- match(betaTrain[,2], rawMAF$SNP) 
              betaTrain$oriSNP <-  rawMAF$SNP[matTrain]
              rownames(betaTrain) <- rawMAF$SNP[matTrain]
              
              
              matTest <- match(betaTest[,2], rawMAF$SNP) 
              betaTest$oriSNP <-  rawMAF$SNP[matTest]
              rownames(betaTest) <- rawMAF$SNP[matTest]
              
              

              ## ================
              ## Accumulate per-CHR
              ## ================
              if (chr == 1) {
                AllBetaTrain      <- betaTrain[,6,drop=FALSE]
                AllBetaTest       <- betaTest[,6,drop=FALSE]
                AllrefAlleleTrain <- betaTrain[,c(7,4)]
                AllrefAlleleTest  <- betaTest[,c(7,4)]
              } else {
                AllBetaTrain      <- rbind(AllBetaTrain, betaTrain[,6,drop=FALSE])
                AllBetaTest       <- rbind(AllBetaTest,  betaTest[,6,drop=FALSE])
                AllrefAlleleTrain <- rbind(AllrefAlleleTrain, betaTrain[,c(7,4)])
                AllrefAlleleTest  <- rbind(AllrefAlleleTest,  betaTest[,c(7,4)])
              }

              rm(betaTrain, betaTest)
            }

            ## Save allele reference files
            refAllelefile_list = list()
            refAllelefile_list[[1]] = paste0(output_sub_folder_predi,"refAllelefile1.txt")
            refAllelefile_list[[2]] = paste0(output_sub_folder_predi,"refAllelefile2.txt")

            write.table(AllrefAlleleTrain,file=refAllelefile_list[[1]],quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
            write.table(AllrefAlleleTest,file=refAllelefile_list[[2]],quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)

            ## Beta matrices
            AllBetaMatrix_list <- list(
              t(AllBetaTrain),
              t(AllBetaTest)
            )

            ## Run R2 + scoring
            GenPreR2_Chrs_saveScore(
              AllBetaMatrix_list,
              savef2_list,
              refAllelefile_list,
              output_sub_folder_predi,
              PV,
              popuseY,
              dirSimuDataSet,
              dirBase,
              sim,
              savename,
              Scoreoutput,
              savecore
            )
          }
        }
      }
    }
  }
}
