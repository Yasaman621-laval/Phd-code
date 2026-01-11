rm(list = ls())

#------------------------------------------------------------
# 0. Helper(s)
#------------------------------------------------------------
count_nonNA <- function(xx) {
  length(which(!is.na(xx)))
}

#------------------------------------------------------------
# 1. Base setup
#------------------------------------------------------------
trait_type <- ""
model      <- "linear"

nsim      <- 10
trait_ids <- c(1, 2, 3)
pop       <- "EUR"

dirBase        <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"

dirSimuTrait   <- paste0(dirBase, "/three_traits/")

#------------------------------------------------------------
# 2. Scenario grid (must match Stage-1 / Stage-4)
#------------------------------------------------------------
#rho_grid <- c(0.2, 0.4, 0.6, 0.8, 0.9)
#hsq_grid <- c(0.3, 0.5, 0.6, 0.8, 0.9)

#rho_grid <- c( 0.4, 0.6, 0.8, 0.9)
#hsq_grid <- c(0.5, 0.6, 0.8, 0.9)

rho_grid = 0.2
hsq_grid = 0.3

#train_grid <- list(50000, 80000, 100000, 150000, 200000)
#val_grid   <- list(5000,  10000, 15000, 20000, 25000)


#train_grid <- list(80000, 100000, 150000, 200000)
#val_grid   <- list(10000, 15000, 20000, 25000)

train_grid = 50000
val_grid = 5000
#------------------------------------------------------------
# 3. Main scenario loop
#------------------------------------------------------------
for (i in seq_along(rho_grid)) {

  rho          <- rho_grid[i]
  hsq          <- hsq_grid[i]
  TrainingNsam <- train_grid[[i]]

  scenario_dir <- paste0(
    dirSimuTrait,
    "sim_hsq", hsq,
    "_rho", rho,
    "_trainEUR", TrainingNsam, "/"
  )

  if (!dir.exists(scenario_dir)) {
    stop("Scenario folder missing: ", scenario_dir,
         "\nRun Stage-1/2/4 first for this scenario.")
  }

  cat("\n================================================\n")
  cat("Stage-5: Build 3-trait GWAS matrices\n")
  cat("Scenario", i,
      "| hsq =", hsq,
      "| rho =", rho,
      "| Train_EUR =", TrainingNsam, "\n")
  cat("Directory:", scenario_dir, "\n")
  cat("================================================\n")

  #--------------------------------------------------------
  # 4. Loop over simulations
  #--------------------------------------------------------
  for (sim in 1:nsim) {

    cat("\n  ---- Simulation", sim, "----\n")

    otherList <- list()
    keepsnps  <- NULL

    #------------------------------------------------------
    # 4.1 Read GWAS results for each trait
    #------------------------------------------------------
    for (iiIndex in trait_ids) {

      N <- TrainingNsam

      fileGWASbeta <- paste0(
        scenario_dir,
        pop, "_trait", iiIndex,
        "_AllChrs_bedformat_sim", sim
      )

      gfile <- paste0(fileGWASbeta, ".assoc.", model)
      if (!file.exists(gfile)) {
        stop("Missing GWAS file: ", gfile,
             "\nRun Stage-4 first.")
      }

      GWASbeta <- read.table(gfile, header = TRUE, as.is = TRUE)
      GWASbeta <- GWASbeta[, c("SNP", "A1", "STAT", "P")]
      names(GWASbeta) <- c("SNP", "A1", "Zobs", "p")

      # Convert Z to beta and SE
      GWASbeta$b  <- GWASbeta$Zobs / sqrt(N)
      GWASbeta$SE <- GWASbeta$b / GWASbeta$Zobs

      # For the first trait, define the base SNP list
      if (iiIndex == trait_ids[1]) {
        keepsnps <- GWASbeta$SNP
      }

      # Suffix all columns by trait index
      names(GWASbeta) <- paste0(names(GWASbeta), iiIndex)
      otherList[[iiIndex]] <- GWASbeta
      rm(GWASbeta)
    }

    nset <- length(trait_ids)

    #------------------------------------------------------
    # 4.2 Merge traits by SNP
    #------------------------------------------------------
    Start <- 1
    for (ii in trait_ids) {

      temp     <- otherList[[ii]]
      nametemp <- paste0("SNP", ii)

      mat <- match(keepsnps, temp[, nametemp])

      if (Start == 1) {
        mainGWAS <- temp[mat, , drop = FALSE]
        Start <- 0
      } else {
        mainGWAS <- cbind(mainGWAS, temp[mat, , drop = FALSE])
      }
    }
    rm(otherList)

    unames <- paste0("Zobs", trait_ids)

    # Quality check (not used, but useful for diagnostics)
    apply(mainGWAS[, unames], 2, summary)

    countsNonNA <- apply(mainGWAS[, unames], 1, count_nonNA)
    table(countsNonNA)

    wkeep <- which(countsNonNA == nset)
    mainGWAS <- mainGWAS[wkeep, ]

    apply(mainGWAS[, paste0("b",  trait_ids)], 2, summary)
    apply(mainGWAS[, paste0("SE", trait_ids)], 2, summary)

    #------------------------------------------------------
    # 4.3 Attach CHR / base SNP / A1
    #------------------------------------------------------
    mainindex <- trait_ids[1]

    mainGWAS$SNP <- mainGWAS[, paste0("SNP", mainindex)]
    mainGWAS$A1  <- mainGWAS[, paste0("A1",  mainindex)]

    # Re-read one GWAS file to get CHR info
    fileGWASbeta0 <- paste0(
      scenario_dir,
      pop, "_trait", mainindex,
      "_AllChrs_bedformat_sim", sim
    )
    snpinfo <- read.table(paste0(fileGWASbeta0, ".assoc.", model),
                          header = TRUE, as.is = TRUE)

    mat <- match(mainGWAS$SNP, snpinfo$SNP)
    if (any(is.na(mat))) {
      cat("Warning: Some SNPs in mainGWAS not found in snpinfo\n")
    }

    snpinfo          <- snpinfo[mat, ]
    mainGWAS$CHR     <- snpinfo$CHR
    mainGWAS$P       <- mainGWAS[, paste0("p", mainindex)]

    #------------------------------------------------------
    # 4.4 Build final 3-trait matrix
    #------------------------------------------------------
    usednames <- c()
    for (ii in trait_ids) {
      usednames <- c(usednames,
                     paste0(c("b", "SE", "p", "Zobs"), ii))
    }

    mainGWAS0 <- mainGWAS[, c("CHR", "SNP", "A1", usednames)]

    cat("  Dimensions of mainGWAS0:", dim(mainGWAS0), "\n")
    if (any(is.na(mainGWAS0))) {
      cat("  Warning: NA values found in mainGWAS0\n")
    }

    #------------------------------------------------------
    # 4.5 Write outputs for this sim + scenario
    #------------------------------------------------------
    out1 <- paste0(
      scenario_dir,
      "GWASbetaStandard_allchrs_3traits_sim", sim, ".txt"
    )
    out2 <- paste0(
      scenario_dir,
      "GWASbeta0Standard_allchrs_3traits_sim", sim, ".txt"
    )

    write.table(mainGWAS0, file = out1,
                quote = FALSE, sep = "\t", row.names = FALSE)
    write.table(mainGWAS,  file = out2,
                quote = FALSE, sep = "\t", row.names = FALSE)

    # Example check: number of SNPs with p2 < 1e-6
    if ("p2" %in% colnames(mainGWAS0)) {
      cat("  # of SNPs with p2 < 1e-6:",
          length(which(mainGWAS0$p2 < 1e-6)), "\n")
    }

    cat("  Columns in mainGWAS0:\n")
    print(colnames(mainGWAS0))

    rm(mainGWAS0, mainGWAS, snpinfo)
    gc()
  } # end sim loop
} # end scenario loop

cat("\n? Stage-5: 3-trait GWAS summary matrices completed for all scenarios and sims.\n")
