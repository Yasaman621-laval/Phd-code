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

nsim = 10
phenoIndex = 3
functionsfolder = "/lustre09/project/6005709/yatah3/Rfunction_PRS_CSX/"
source(paste0(functionsfolder,"PRS_utility.r"))
source(paste0(functionsfolder,"Iterative_Rfunctions.r"))
source(paste0(functionsfolder,"AllPRS_Rfunctions.r"))

dirtemp = Sys.getenv('SLURM_TMPDIR')
folder_predi = paste0(dirtemp,"/")

R2output = paste0(sc_dir,"/",savename,"/R2output_meta/")
system(paste0("mkdir -p ",R2output))
Scoreoutput = paste0(folder_predi,savename,"/Scoreoutput_meta/")
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

for (sim in 6:nsim) {

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

      files_meta    <- finishedfiles_all[ grepl("_META_", finishedfiles_all) ]
      

      ############################################################
      ##  B) META BRANCH (1 file per chromosome)
      ############################################################
      if (length(files_meta) > 0) {

        finishedchr_meta <- c()
        for (fn in files_meta) {
          parts   <- strsplit(fn, "_")[[1]]
          chrpart <- parts[grep("chr", parts)]
          chrnum  <- as.numeric(gsub("chr|\\.txt","", chrpart))
          finishedchr_meta <- c(finishedchr_meta, chrnum)
        }

        test_meta <- table(finishedchr_meta)

        ## Expect exactly 1 META output per chr
        if (length(which(test_meta == 1)) == 22) {

          savecore_meta <- paste0(popTrain, "_", pop,
                                  "sim", sim, "phiuse", phiuse, "_META")

          savef2_list_meta <- list(
            paste0(R2output, savecore_meta, "_Meta_results.RData")
          )

          AllFinish_meta = 1
          if (file.exists(savef2_list_meta[[1]])) {
            load(savef2_list_meta[[1]])
            if (exists("Finish") && Finish == 0) AllFinish_meta = 0
          } else AllFinish_meta = 0

          if (AllFinish_meta == 0) {

            output_sub_folder_predi_meta <-
              paste0(folder_predi, savecore_meta, "/")
            system(paste0("mkdir -p ", output_sub_folder_predi_meta))

           for (chr in 1:22) {

              chrname = paste0("chr", chr, ".txt")
              idx     = grep(chrname, files_meta)
              usefiles = files_meta[idx]

              if (length(usefiles) != 1)
                stop("Expected 1 NON-META files for chr=", chr,
                     " but found ", length(usefiles))

              extract_pop <- function(filename) {
                # captures pop between _simX_ and _pst_eff
                sub(".*_sim[0-9]+_([^_]+)_pst_eff.*", "\\1", filename)
              }

           
              file = paste0(output_sub_folder, usefiles)
         
              beta <- read.table(file, as.is=TRUE)
          

              ## ================
              ## Correct mapping
              ## ================
              mat <- match(beta[,2], rawMAF$SNP) 
              beta$oriSNP <-  rawMAF$SNP[mat]
              rownames(beta) <- rawMAF$SNP[mat]
              
              

              ## ================
              ## Accumulate per-CHR
              ## ================
              if (chr == 1) {
                AllBeta      <- beta[,6,drop=FALSE]
                AllrefAllele <- beta[,c(7,4)]
              } else {
                AllBeta    <- rbind(AllBeta, beta[,6,drop=FALSE])
                AllrefAllele <- rbind(AllrefAllele, beta[,c(7,4)])
              }

              rm(beta)
            }

            ## Save allele reference files
            refAllelefile_list = list()
            refAllelefile_list[[1]] = paste0(output_sub_folder_predi_meta,"refAllelefile.txt")

            write.table(AllrefAllele,file=refAllelefile_list[[1]],quote = FALSE, sep = "\t",row.names = FALSE,col.names = FALSE)
       

            ## Beta matrices
            AllBetaMatrix_list <- list(
              t(AllBeta)
            )

            ## Run R2 + scoring
            GenPreR2_Chrs_saveScore(
              AllBetaMatrix_list,
              savef2_list_meta,
              refAllelefile_list,
              output_sub_folder_predi_meta,
              PV,
              popuseY,
              dirSimuDataSet,
              dirBase,
              sim,
              savename,
              Scoreoutput,
              savecore_meta
            )
          }
        }
      }
    }
  }
}
