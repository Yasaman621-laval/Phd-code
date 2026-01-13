rm(list=ls())
library(data.table)
library(dplyr)

#------------------------------------------------------------#
# DIRECTORIES – EDIT THIS ONLY IF NEEDED
#------------------------------------------------------------#
dirBase   <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype"
out_file  <- "/lustre06/project/6005709/yatah3/simulation/project2/PRs-csx/snpinfo_mult_1kg_hm3"

popvec  <- c("EUR","EAS","AFR")

#------------------------------------------------------------#
# Function to read BIM + allele frequencies for POP
#------------------------------------------------------------#
read_pop_info <- function(pop) {
  
  prefix <- paste0(dirBase, "/", pop, "AllChrs_bedformat")
  bim_file <- paste0(prefix, ".bim")
  
  bim <- fread(bim_file, header = FALSE)
  colnames(bim) <- c("CHR","SNP","CM","BP","A1","A2")
  
  # Compute allele frequencies using PLINK
  freq_prefix <- tempfile()
  cmd <- paste("plink --bfile", prefix, "--freq --out", freq_prefix)
  system(cmd)
  
  frq <- fread(paste0(freq_prefix, ".frq"))
  colnames(frq) <- c("CHR","SNP","A1_frq","A2_frq","AF","NCHROBS")
  
  merged <- merge(bim, frq[,c("SNP","AF")], by="SNP", all.x=TRUE)
  merged <- merged %>% arrange(CHR, BP)
  
  return(merged)
}

#------------------------------------------------------------#
# Read EUR, EAS, AFR
#------------------------------------------------------------#
eur <- read_pop_info("EUR")
eas <- read_pop_info("EAS")
afr <- read_pop_info("AFR")

#------------------------------------------------------------#
# Keep only SNPs present in ALL populations
#------------------------------------------------------------#
common_snps <- Reduce(intersect, list(eur$SNP, eas$SNP, afr$SNP))

eur <- eur[eur$SNP %in% common_snps, ]
eas <- eas[eas$SNP %in% common_snps, ]
afr <- afr[afr$SNP %in% common_snps, ]

#------------------------------------------------------------#
# Build unified SNPINFO (use EUR for CHR/BP/A1/A2)
#------------------------------------------------------------#
df <- data.frame(
  CHR = eur$CHR,
  SNP = eur$SNP,
  BP  = eur$BP,
  A1  = eur$A1,
  A2  = eur$A2,
  
  FRQ_EUR = eur$AF,
  FRQ_EAS = eas$AF,
  FRQ_AFR = afr$AF
)

#------------------------------------------------------------#
# Create PRS-CSx FULL 15-column SNPINFO
#------------------------------------------------------------#
df_out <- data.frame(
  CHR = df$CHR,
  SNP = df$SNP,
  BP  = df$BP,
  A1  = df$A1,
  A2  = df$A2,

  FRQ_AFR = df$FRQ_AFR,
  FRQ_AMR = NA,             # placeholder
  FRQ_EAS = df$FRQ_EAS,
  FRQ_EUR = df$FRQ_EUR,
  FRQ_SAS = NA,             # placeholder

  FLP_AFR = ifelse(is.na(df$FRQ_AFR), 0, 1),
  FLP_AMR = 0,
  FLP_EAS = ifelse(is.na(df$FRQ_EAS), 0, 1),
  FLP_EUR = ifelse(is.na(df$FRQ_EUR), 0, 1),
  FLP_SAS = 0
)

#------------------------------------------------------------#
# Save SNP info file
#------------------------------------------------------------#
fwrite(df_out, out_file, sep="\t", quote=FALSE, row.names=FALSE)
cat("Created:", out_file, "with", nrow(df_out), "SNPs\n")
