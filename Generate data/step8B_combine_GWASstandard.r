rm(list=ls())
#module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
#module load intel/2016.4 r/3.5.0

trait_type = ""
model = "linear"

count_nonNA = function(xx){
  length(which(!is.na(xx)))
}

iiset = c(1,2,3)
nsim = 10
popvec = c("AFR","EAS","EUR")

#-------------------------------------------------------------
# Base directory holding all scenario folders
#-------------------------------------------------------------
dirBase = "/lustre06/project/6005709/yatah3/simulation/SimuGenotype"
scenario_dirs = list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs = scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs)==0) stop("? No scenario folders found under ", dirBase)

#-------------------------------------------------------------
# Loop through each scenario
#-------------------------------------------------------------
for (sc_dir in scenario_dirs) {

  sc_name = basename(sc_dir)
  cat("\n=====================================\n")
  cat("Processing scenario:", sc_name, "\n")
  cat("=====================================\n")

  #-----------------------------------------------
  # Parse training sample sizes from folder name
  #-----------------------------------------------
  train_str = sub(".*_train", "", sc_name)
  train_str_clean = gsub("e\\+0", "e", train_str)
  nums = as.numeric(sapply(strsplit(train_str_clean, "-")[[1]], function(x) eval(parse(text = x))))
  TrainingNsam = nums

  cat("TrainingNsam =", paste(TrainingNsam, collapse=", "), "\n")

  dirSimuData = paste0(sc_dir, "/")

  #-----------------------------------------------
  # Loop through each simulation
  #-----------------------------------------------
  for (sim in 1:nsim) {
    print(paste("Simulation:", sim))
    otherList = list()

    #---------------------------------------------
    # Step 1. Read GWAS for each population
    #---------------------------------------------
    for (iiIndex in iiset) {
      N = TrainingNsam[iiIndex]
      pop = popvec[iiIndex]
      fileGWASbeta = paste0(dirSimuData, pop, "_assoc_linear_", "sim", sim)
      infile = paste0(fileGWASbeta, ".assoc.", model)

      if (!file.exists(infile)) {
        warning("Missing file: ", infile)
        next
      }

     GWASbeta <- read.table(
  infile,
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "",      # any number of spaces counts as a separator
  fill = TRUE,   # pad missing columns with NA instead of failing
  strip.white = TRUE # remove extra spaces around values
)

      GWASbeta = GWASbeta[,c("SNP","A1","STAT","P")]
      names(GWASbeta) = c("SNP","A1","Zobs","p")
      GWASbeta$b = GWASbeta$Zobs / sqrt(N)
      GWASbeta$SE = GWASbeta$b / GWASbeta$Zobs

      if (iiIndex==1) {
        keepsnps = GWASbeta$SNP
      }

      names(GWASbeta) = paste0(names(GWASbeta), iiIndex)
      otherList[[iiIndex]] = GWASbeta
      rm(GWASbeta)
    }

    nset = length(iiset)

    #---------------------------------------------
    # Step 2. Merge across populations
    #---------------------------------------------
    Start = 1
    for (ii in iiset) {
      temp = otherList[[ii]]
      nametemp = paste0("SNP", ii)
      mat = match(keepsnps, temp[, nametemp])
      if (Start==1) {
        mainGWAS = temp[mat,,drop=FALSE]
        Start = 0
      } else {
        mainGWAS = cbind(mainGWAS, temp[mat,,drop=FALSE])
      }
    }
    rm(otherList)

    unames = paste0("Zobs", iiset)
    countsNonNA = apply(mainGWAS[,unames], 1, count_nonNA)
    wkeep = which(countsNonNA==nset)
    mainGWAS = mainGWAS[wkeep,]

    #---------------------------------------------
    # Step 3. Attach CHR and P info
    #---------------------------------------------
    mainindex = 1
    mainGWAS$SNP = mainGWAS[, paste0("SNP", mainindex)]
    mainGWAS$A1  = mainGWAS[, paste0("A1", mainindex)]

    pop0 = popvec[mainindex]
    fileGWASbeta0 = paste0(dirSimuData, pop0, "_assoc_linear_", "sim", sim)
    snpinfo <- read.table(
  file = paste0(fileGWASbeta0, ".assoc.", model),
  header = TRUE,
  stringsAsFactors = FALSE,
  sep = "",       # treat any number of spaces as separator
  fill = TRUE,    # pad missing values with NA instead of throwing error
  strip.white = TRUE # trim extra spaces
)


    mat = match(mainGWAS$SNP, snpinfo$SNP)
    snpinfo = snpinfo[mat,]
    mainGWAS$CHR = snpinfo$CHR
    mainGWAS$P = mainGWAS[, paste0("p", mainindex)]

    #---------------------------------------------
    # Step 4. Select final columns
    #---------------------------------------------
    usednames = c()
    for (ii in iiset) {
      usednames = c(usednames, paste0(c("b","SE","p","Zobs"), ii))
    }

    mainGWAS0 = mainGWAS[, c("CHR","SNP","A1", usednames)]

    #---------------------------------------------
    # Step 5. Save output (overwrite existing)
    #---------------------------------------------
    out1 = paste0(dirSimuData, "GWASbetaStandard_allchrs", "sim", sim, ".txt")
    out2 = paste0(dirSimuData, "GWASbeta0Standard_allchrs", "sim", sim, ".txt")

    if (file.exists(out1)) file.remove(out1)
    if (file.exists(out2)) file.remove(out2)

    write.table(mainGWAS0, file=out1, quote=FALSE, sep="\t", row.names=FALSE)
    write.table(mainGWAS,  file=out2, quote=FALSE, sep="\t", row.names=FALSE)

    print(paste("Saved:", out1))
    rm(mainGWAS0, mainGWAS, snpinfo)
    gc()
  }

  cat("Completed scenario:", sc_name, "\n")
}

cat("\n? All scenarios processed successfully.\n")
