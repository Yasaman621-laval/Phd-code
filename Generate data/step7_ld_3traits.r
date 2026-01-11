rm(list = ls())
library(MASS)
library(mvtnorm)
library(gtools)

#------------------------------------------------------------
# 1. Base setup
#------------------------------------------------------------
popvec <- c("AFR", "EAS", "EUR")
trait_ids <- 1:3

dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"

dirSimuTraitOut <- paste0(dirBase, "/three_traits/")

# SNP list for extraction
totalsfile <- paste0(dirBase, "SNPs_allChrs.txt")
AllSNPs <- read.table(totalsfile, as.is = TRUE)[, 1]
usedsnpfile <- totalsfile

#------------------------------------------------------------
# 2. Reference allele file (only created once)
#------------------------------------------------------------
refAllelefile <- paste0(dirBase, "refereceAlleleSNPs.txt")

if (!file.exists(refAllelefile)) {
  pop_ref <- "EUR"   # use EUR as reference
  plink_bim <- paste0(dirSimuTraitOut, pop_ref, "AllChrs_bedformat.bim")

  if (!file.exists(plink_bim)) {
    stop("Reference BIM file missing: ", plink_bim)
  }

  bim <- read.table(plink_bim, as.is = TRUE)
  refSNP <- bim[, c(2, 5)]  # SNP ID + reference allele
  write.table(refSNP, file = refAllelefile, quote = FALSE, sep = "\t",
              row.names = FALSE, col.names = FALSE)

  rm(bim)
}

cat("Reference allele file:", refAllelefile, "\n")

#------------------------------------------------------------
# 3. LD computation parameters
#------------------------------------------------------------
numsnp <- 400     # LD window size
pop <- "EUR"       # Only EUR LD window needed (as in your code)

#------------------------------------------------------------
# 4. LD Calculation Loop
#    populations × traits × chromosomes
#------------------------------------------------------------
for (pop in popvec) {

  cat("\n=============================================\n")
  cat("Computing LD for population:", pop, "\n")
  cat("=============================================\n")

  for (trait in trait_ids) {

    cat("\n --- Trait", trait, "---\n")

    for (chr in 1:22) {

      cat("      Chromosome", chr, "\n")

      # Folder for the chromosome
      chrDir <- paste0(dirSimuTraitOut, "Chr", chr, "/")
      dirPlinkFormat <- paste0(chrDir, "PlinkFormat/")
      system(paste0("mkdir -p ", dirPlinkFormat))

      # Original PLINK file for this chr/pop
      plinkfile <- paste0(dirPlinkFormat, pop, "chr", chr, "bedformat")

      if (!file.exists(paste0(plinkfile, ".bim"))) {
        stop("Missing PLINK chr file: ", plinkfile)
      }

      # Output name
      outname <- paste0(
        dirPlinkFormat,
        pop, "_trait", trait, "_chr", chr,
        "_bedformat_numsnp", numsnp
      )

      # LD command
      cmd <- paste0(
        "plink --noweb --allow-no-sex",
        " --bfile ", plinkfile,
        " --reference-allele ", refAllelefile,
        " --extract ", usedsnpfile,
        " --r",
        " --ld-window ", numsnp,
        " --out ", outname
      )

      system(cmd)
    }
  }
}

cat("\n? Stage-3 LD computation complete for all populations, traits and 22 chromosomes.\n")
