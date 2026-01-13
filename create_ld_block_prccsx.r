#!/usr/bin/env Rscript
rm(list=ls())

library(data.table)
library(hdf5r)

cat("\n===== MEMORY-SAFE PRS-CSx LD BLOCK BUILDER =====\n")

# ------------------------------------------------------------
# Argument
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript ld_builder_stream.R <jset>")
jset <- as.numeric(args[1])

# ------------------------------------------------------------
# Input configuration
# ------------------------------------------------------------
input_file <- "/lustre06/project/6005709/yatah3/simulation/project2/PRs-csx/input_ld.txt"
inputs <- read.table(input_file, header=TRUE, as.is=TRUE)

chrom    <- inputs$chr[jset]
popIndex <- inputs$pop[jset]
popvec   <- c("AFR","EAS","EUR")
pop      <- popvec[popIndex]

cat("Chromosome:", chrom, "Population:", pop, "\n")

# ------------------------------------------------------------
# LD file + BIM
# ------------------------------------------------------------
ld_file <- sprintf(
  "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/sim_hsq0.5_rho0.4_train10000-10000-80000/Chr%d/PlinkFormat/%schr%dbedformat_numsnp400.ld",
  chrom, pop, chrom
)

bim_file <- sprintf(
  "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/%sAllChrs_bedformat.bim",
  pop
)

# ------------------------------------------------------------
# Read inputs
# ------------------------------------------------------------
bim <- fread(bim_file, header=FALSE)
colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
bim <- bim[bim$CHR == chrom]
bim <- bim[order(BP)]

snps <- bim$SNP
positions <- bim$BP
n <- length(snps)

cat("Total SNPs:", n, "\n")

# ------------------------------------------------------------
# Sliding window block builder (1 Mb)
# ------------------------------------------------------------
block_list <- list()
block_idx <- 1
start <- 1

while (start <= n) {
  end <- start
  while (end < n && positions[end+1] - positions[start] <= 1e6)
    end <- end + 1
  block_list[[block_idx]] <- c(start, end)
  block_idx <- block_idx + 1
  start <- end + 1
}

cat("Total blocks:", length(block_list), "\n")

# ------------------------------------------------------------
# Output directory
# ------------------------------------------------------------
out_dir <- sprintf(
  "/lustre06/project/6005709/yatah3/simulation/project2/PRs-csx/ldblk_1kg_%s",
  tolower(pop)
)

dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

outfile <- sprintf("%s/ldblk_1kg_chr%d.hdf5", out_dir, chrom)

if (file.exists(outfile))
  file.remove(outfile)

# ------------------------------------------------------------
# Prepare LD
# ------------------------------------------------------------
ld <- fread(ld_file)

ld[, i := match(SNP_A, snps)]
ld[, j := match(SNP_B, snps)]
ld <- ld[!is.na(i) & !is.na(j)]

# ------------------------------------------------------------
# Create HDF5
# ------------------------------------------------------------
h5 <- H5File$new(outfile, mode="w")

for (b in seq_along(block_list)) {
  idx <- block_list[[b]]
  s <- idx[1]; e <- idx[2]
  k <- e - s + 1

  cat("Block", b, " size =", k, "\n")

  LD_block <- diag(1, k)

  ld_sub <- ld[i >= s & i <= e & j >= s & j <= e]
  ld_sub[, ii := i - s + 1]
  ld_sub[, jj := j - s + 1]

  for (row in 1:nrow(ld_sub)) {
    ii <- ld_sub$ii[row]
    jj <- ld_sub$jj[row]
    r  <- ld_sub$R[row]
    LD_block[ii, jj] <- r
    LD_block[jj, ii] <- r
  }

  grp <- h5$create_group(paste0("blk_", b))
  
grp$create_dataset(
  name = "ldblk",
  robj = LD_block
)

grp$create_dataset(
  name = "snplist",
  robj = as.character(snps[s:e])
)

  grp$close()
}

h5$close_all()

cat("\n===== DONE: MEMORY-SAFE PRS-CSx LD BUILDER =====\n")
