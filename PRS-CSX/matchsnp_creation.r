############################################################
# Generate matchSNP reference file for PRS-CSx
# Input: BIM + FRQ files (AFR/EUR/EAS)
# Output: matchSNP_<POP>.txt for each population
############################################################
rm(list=ls())
library(dplyr)
library(readr)

# ======= 1. Set base directory =======
dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"

setwd(dirBase)

# ======= 2. Detect BIM and FRQ files =======
bim_files <- list.files(pattern = "_bedformat.bim$")
frq_files <- list.files(pattern = "_bedformat.frq$")

cat("Found BIM files:\n")
print(bim_files)

cat("Found FRQ files:\n")
print(frq_files)

# ======= 3. Function to process one population =======
make_match_snp <- function(pop){

  cat("\n==============================\n")
  cat("Processing population:", pop, "\n")
  cat("==============================\n")

  bim_file <- paste0(pop, "AllChrs_bedformat.bim")
  frq_file <- paste0(pop, "AllChrs_bedformat.frq")

  if (!file.exists(bim_file)) stop("Missing BIM file: ", bim_file)
  if (!file.exists(frq_file)) stop("Missing FRQ file: ", frq_file)

  # --- 3a. Load BIM ---
  bim <- read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
  # BIM columns:
  # V1 = CHR, V2 = SNP, V3 = CM, V4 = BP, V5 = A1, V6 = A2
  bim <- bim %>% select(CHR = V1, SNP = V2, A1 = V5, A2 = V6)

  # --- 3b. Load FRQ ---
  frq <- read.table(frq_file, header = TRUE, stringsAsFactors = FALSE)
  # Ensure FRQ has: CHR, SNP, A1, A2, MAF, NCHROBS
  frq <- frq %>% select(CHR, SNP, A1, A2, MAF, NCHROBS)

  # --- 3c. Merge BIM + FRQ ---
  merged <- bim %>%
  left_join(frq, by = c("CHR", "SNP", "A1", "A2"))


  if (nrow(merged) == 0){
    stop("No overlapping SNPs between BIM and FRQ for ", pop)
  }

  cat("Merged SNPs:", nrow(merged), "\n")

  # --- 3d. Add OA column ---
  merged <- merged %>%
    mutate(OA = paste0(A1, A2)) %>%
    select(CHR, SNP, A1, A2, MAF, NCHROBS, OA)

  # --- 3e. Output ---
  outFile <- paste0("matchSNP_", pop, ".txt")
  write.table(merged, outFile, quote = FALSE, sep = "\t", row.names = FALSE)

  cat("Output written to:", outFile, "\n")
}

# ======= 4. Run for AFR/EUR/EAS automatically =======
pops <- c("AFR", "EUR", "EAS")

for (p in pops){
  make_match_snp(p)
}

############################################################
# Create final common SNP file across AFR / EUR / EAS
############################################################

library(dplyr)
library(readr)

dirBase <- "/lustre06/project/6005709/yatah3/simulation/SimuGenotype/"
setwd(dirBase)

# 1. Read matchSNP files
afr <- read.table("matchSNP_AFR.txt", header = TRUE, stringsAsFactors = FALSE)
eur <- read.table("matchSNP_EUR.txt", header = TRUE, stringsAsFactors = FALSE)
eas <- read.table("matchSNP_EAS.txt", header = TRUE, stringsAsFactors = FALSE)

# 2. Find common SNPs by ID
common_snps <- Reduce(intersect, list(afr$SNP, eur$SNP, eas$SNP))
cat("Common SNP count:", length(common_snps), "\n")

# 3. Filter each population to common SNPs
afr_c <- afr %>% filter(SNP %in% common_snps)
eur_c <- eur %>% filter(SNP %in% common_snps)
eas_c <- eas %>% filter(SNP %in% common_snps)

# 4. Check allele consistency
check_alleles <- function(df1, df2){
  all(df1$A1 == df2$A1 & df1$A2 == df2$A2)
}

if (!check_alleles(afr_c, eur_c)) stop("Allele mismatch between AFR and EUR")
if (!check_alleles(afr_c, eas_c)) stop("Allele mismatch between AFR and EAS")

# 5. Merge all info into one table
final <- afr_c %>%
  select(CHR, SNP, A1, A2, OA,
         MAF_AFR = MAF, NCHROBS_AFR = NCHROBS) %>%
  left_join(eur_c %>% 
              select(SNP, MAF_EUR = MAF, NCHROBS_EUR = NCHROBS),
            by = "SNP") %>%
  left_join(eas_c %>% 
              select(SNP, MAF_EAS = MAF, NCHROBS_EAS = NCHROBS),
            by = "SNP")

# 6. Order by chromosome and position if SNP is chr:pos
if (grepl(":", final$SNP[1])) {
  final <- final %>%
    mutate(POS = as.numeric(sub(".*:(\\d+)_.*", "\\1", SNP))) %>%
    arrange(CHR, POS) %>%
    select(-POS)
}

# 7. Write output
write.table(final, "matchSNP_allPops.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nFinal merged file written to matchSNP_allPops.txt\n")

