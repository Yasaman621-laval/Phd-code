# =========================
# 0) Setup
# =========================
rm(list = ls())

  library(dplyr)



# ---- PATHS (EDIT) ----
out_dir = "/lustre09/project/6005709/yatah3/real-data/out/"


# PLINK refs (BED/BIM/FAM triple per pop)
plink_dir <- "/lustre09/project/6005709/yatah3/real-data"
bfile_EUR <- file.path(plink_dir, "1kg_hg38_EUR")
bfile_AFR <- file.path(plink_dir, "1kg_hg38_AFR")
bfile_EAS <- file.path(plink_dir, "1kg_hg38_EAS")

plink_bin <- "plink"  # full path if needed

# =========================
# 1) Read gwas from each pop 
# =========================
# === Load EUR data ===
eur <- read.delim("/lustre09/project/6005709/yatah3/real-data/EUR_harmonized.tsv")


eur <- eur %>%
    rename(
    CHR = chromosome,
    POS = base_pair_location,
    Effect_Allele = effect_allele,
    Other_Allele = other_allele,
    beta_EUR = beta_logOR,
    se_EUR = se_recalc,
    Z_EUR= Z_from_pval,
    pval_EUR = p_value,
    EAF_EUR = effect_allele_frequency,
    N_EUR = N
  )%>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele, beta_EUR, se_EUR, pval_EUR, Z_EUR, N_EUR, EAF_EUR)

 



# === Load EAS data ===
eas <- read.delim("/lustre09/project/6005709/yatah3/real-data/EAS_harmonized.tsv")


eas <- eas %>%
  rename(
    CHR = chromosome,
    POS = base_pair_location,
    Effect_Allele = effect_allele,
    Other_Allele = other_allele,
    beta_EAS = beta_logOR,                 # temporary name
    se_EAS = se_recalc,
    pval_EAS = p_value,
    N_EAS = n,
    Z_EAS= Z_from_pval,
     EAF_EAS = effect_allele_frequency
  ) %>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele,beta_EAS, se_EAS, pval_EAS,Z_EAS, N_EAS,EAF_EAS)
 



# === Load AFR data ===
afr <- read.delim("/lustre09/project/6005709/yatah3/real-data/AFR_harmonized.tsv")


# Make sure allele columns are capitalized
afr <- afr %>%
  rename(
    CHR = chromosome,
    POS = base_pair_location,
    Effect_Allele = effect_allele,
    Other_Allele = other_allele,
     beta_AFR = beta_logOR,              # temporary name
    se_AFR = se_recalc,
    pval_AFR = p_value,
    N_AFR = n,
     Z_AFR= Z_from_pval,
     EAF_AFR = effect_allele_frequency) %>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele, beta_AFR, se_AFR, pval_AFR,Z_AFR, N_AFR,EAF_AFR)       
  
  
  
  
 cat("Loaded EUR GWAS rows:", nrow(eur), "\n")
cat("Loaded EAS GWAS rows:", nrow(eas), "\n")
cat("Loaded AFR GWAS rows:", nrow(afr), "\n")
 
         
#___________________________________________________________________

# === Merge all three files on the key ===
# Final merged data
merged <- eur %>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele,
         beta_EUR, se_EUR, pval_EUR,N_EUR,EAF_EUR,Z_EUR) %>%
  inner_join(
    eas %>% select(rsid, beta_EAS, se_EAS, pval_EAS, N_EAS,Z_EAS,EAF_EAS),
    by = "rsid"
  ) %>%
  inner_join(
    afr %>% select(rsid, beta_AFR, se_AFR, pval_AFR, N_AFR,Z_AFR,EAF_AFR),
    by = "rsid"
  )


# === Save to file ===
write.table(
  merged,
  file = file.path(out_dir, "merged_summary_stats.txt"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)
#
cat("SNPs shared across all three populations (pre-QC):", nrow(merged), "\n")



#merged <- read.delim(file.path(gwas_dir, "merged_summary_stats.txt"), sep = "\t", header = TRUE) 

# =========================
# 2) GWAS QC: remove ambiguous + MAF tails across all pops + drop NA pvals
# =========================
ambig_pairs <- c("AT","TA","CG","GC")
merged <- merged %>%
  filter(!(paste0(Effect_Allele, Other_Allele) %in% ambig_pairs))

# MAF (effect allele freq) between 0.05 and 0.95 across pops
merged <- merged %>%
  filter(
    !is.na(EAF_EUR) & !is.na(EAF_EAS) & !is.na(EAF_AFR),
    EAF_EUR > 0.05 & EAF_EUR < 0.95,
    EAF_EAS > 0.05 & EAF_EAS < 0.95,
    EAF_AFR > 0.05 & EAF_AFR < 0.95
  )

# Remove NA p-values
merged <- merged %>%
  filter(!is.na(pval_EUR), !is.na(pval_EAS), !is.na(pval_AFR))

# Optional: ensure numeric
merged <- merged %>%
  mutate(
    CHR = as.integer(CHR),
    POS = as.integer(POS)
  )

 write.table(merged, file.path(out_dir, "merged_summary_stats.filtered.tsv"), sep = "\t")



cat("SNPs remaining after QC and MAF filters:", nrow(merged), "\n")


# =========================
# 3) Build a shared SNP backbone from BIMs (EUR/AFR/EAS)
#    - keep biallelic SNVs (allele length == 1)
#    - drop ambiguous AT/TA/CG/GC
#    - intersect across EUR/AFR/EAS
#    - intersect with GWAS rsids
# =========================
read_bim_clean <- function(bim_path) {
  bim <- read.delim(paste0(bim_path, ".bim"), header = FALSE)
 colnames(bim) <- c("CHR", "SNP", "CM", "POS", "A1", "A2")
  # keep SNVs only
 bim <- bim[nchar(bim$A1) == 1 & nchar(bim$A2) == 1, ]
  amb <- paste0(toupper(bim$A1), toupper(bim$A2))
bim <- bim[!(amb %in% ambig_pairs), ]

}

bEUR <- read_bim_clean(bfile_EUR)
bAFR <- read_bim_clean(bfile_AFR)
bEAS <- read_bim_clean(bfile_EAS)

# rsid intersection across reference panels
backbone_rsids <- Reduce(intersect, list(bEUR$SNP, bAFR$SNP, bEAS$SNP))

# Also require presence in your filtered merged GWAS
#based on the seperate gwas,r , we got this is the right order for allele when copare with bim file:
merged <- merged %>% mutate(SNPID = paste0(CHR, ":", POS, ":", Other_Allele, ":", Effect_Allele))

backbone_rsids <- intersect(backbone_rsids, merged$SNPID)
length(backbone_rsids)

# Write full backbone (all chromosomes mixed; PLINK can still use it with --extract)
backbone_file <- file.path(out_dir, "SNPs_allChrs.txt")
write.table(backbone_rsids,
            file = backbone_file,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
cat("Final backbone SNPs (shared across reference panels):", length(backbone_rsids), "\n")


# =========================
# 4) Make per-chromosome SNP–P files (per pop) restricted to backbone
#    (PLINK --clump needs 'SNP    P')
# =========================
mk_snp_p_files <- function(df, pcol, pop_label) {
  df2 <- df %>%
    filter(SNPID %in% backbone_rsids) %>%
    select(SNPID, CHR, P = !!sym(pcol)) %>%
    distinct(SNPID, .keep_all = TRUE)

  for (chr in 1:22) {
    tmp <- df2 %>% filter(CHR == chr) %>% select(SNP = SNPID, P)
    if (nrow(tmp) == 0) next
    outf <- file.path(out_dir, sprintf("SNPs_pval_chr%d_%s.txt", chr, pop_label))
            
     write.table(tmp,
              file = outf,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)         
            
  }
}

mk_snp_p_files(merged, "pval_EUR", "EUR")
mk_snp_p_files(merged, "pval_AFR", "AFR")
mk_snp_p_files(merged, "pval_EAS", "EAS")

# =========================
# 5) reference-allele files  (SNPID, A1)
# =========================
outpath = file.path(out_dir, "reference_alleles.txt")
ref_file = merged %>%   
          filter(SNPID %in% backbone_rsids) %>%
    transmute(SNPID = SNPID, A1 = Effect_Allele) %>%
    distinct(SNPID, .keep_all = TRUE)

write.table(ref_file, outpath, sep = "\t", col.names = FALSE)


# =========================
# 6) Run PLINK clumping per pop on the same backbone
#    (Do NOT intersect post-clump lead SNPs across pops)
# =========================
r2cut=0.05

run_clump <- function(bfile, pop_label) {
  for (chr in 1:22) {
    snp_p <- file.path(out_dir, sprintf("SNPs_pval_chr%d_%s.txt", chr, pop_label))
    if (!file.exists(snp_p)) next
    outpref <- file.path(out_dir, sprintf("clumped_chr%d_r2cut0.05_%s", chr, pop_label))

    
     command = paste0("plink --noweb --allow-no-sex --bfile ", bfile , " --extract ", backbone_file, " --clump ", snp_p," --clump-p1 1.0 --clump-p2 1.0 --clump-r2 ",r2cut," --clump-kb 500",
	     " --out ", outpref);

system(command)  
    
  }
}

run_clump(bfile_EUR, "EUR")
run_clump(bfile_AFR, "AFR")
run_clump(bfile_EAS, "EAS")

# =========================
# 7) Collect clumped index SNPs per pop and export "final" clumped GWAS per pop
# =========================
collect_clumped_snps <- function(pop_label) {
  kept <- c()
  for (chr in 1:22) {
    f <- file.path(out_dir, sprintf("clumped_chr%d_r2cut0.05_%s.clumped", chr, pop_label))
    if (!file.exists(f)) next

    # Safely try reading the file
    d <- tryCatch(read.table(f, header = TRUE, stringsAsFactors = FALSE),
                  error = function(e) NULL)
    if (is.null(d) || nrow(d) == 0) next

    # PLINK marks index SNPs with F == 1
    if ("F" %in% colnames(d)) {
      d <- d[d$F == 1, , drop = FALSE]
      if (nrow(d) == 0) next

      # save per-chromosome SNP list (needed for LD)
      subsnps <- d$SNP
      write.table(
        subsnps,
        file = file.path(out_dir, paste0("Chr", chr, "r2cut", 0.05, pop_label, ".txt")),
        quote = FALSE,
        sep = "\t",
        row.names = FALSE,
        col.names = FALSE
      )

      kept <- c(kept, subsnps)
    } else {
      warning(sprintf("Column 'F' not found in %s", basename(f)))
    }
  }
  unique(kept)
}



kept_EUR <- collect_clumped_snps("EUR")
kept_AFR <- collect_clumped_snps("AFR")
kept_EAS <- collect_clumped_snps("EAS")

cat("EUR lead SNPs after clumping:", length(kept_EUR), "\n")
cat("AFR lead SNPs after clumping:", length(kept_AFR), "\n")
cat("EAS lead SNPs after clumping:", length(kept_EAS), "\n")
cat("Union of all lead SNPs across populations:", length(all_kept_snps), "\n")


# Final per-pop clumped GWAS tables (NO intersection)
final_EUR <- merged %>% filter(SNPID %in% kept_EUR) %>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele,
         beta_EUR, se_EUR, pval_EUR, Z_EUR, N_EUR, EAF_EUR,SNPID)
final_AFR <- merged %>% filter(SNPID %in% kept_AFR) %>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele,
         beta_AFR, se_AFR, pval_AFR, Z_AFR, N_AFR, EAF_AFR,SNPID)
final_EAS <- merged %>% filter(SNPID %in% kept_EAS) %>%
  select(rsid, CHR, POS, Effect_Allele, Other_Allele,
         beta_EAS, se_EAS, pval_EAS, Z_EAS, N_EAS, EAF_EAS,SNPID)

# EUR
write.table(final_EUR,
            file = file.path(out_dir, "final_clumped_EUR.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# AFR
write.table(final_AFR,
            file = file.path(out_dir, "final_clumped_AFR.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# EAS
write.table(final_EAS,
            file = file.path(out_dir, "final_clumped_EAS.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# ---------- union of SNPs across pops ----------
all_kept_snps <- unique(c(final_EUR$SNPID, final_AFR$SNPID, final_EAS$SNPID))
length(all_kept_snps)



# ---------- subset merged GWAS to that union ----------
final_cross <- merged %>% filter(SNPID %in% all_kept_snps)

# ---------- order by chromosome & position ----------
final_cross <- final_cross %>% arrange(as.numeric(CHR), as.numeric(POS))

# ---------- save ----------
outfile <- file.path(out_dir, "final_clumped_crosspop.tsv")
write.table(final_cross,
            file = outfile,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


cat("✅ Cross-ancestry final GWAS saved as:", outfile, "\n")
cat("Total SNPs:", nrow(final_cross), "\n")



# =====================================================
# Step 8: compute LD
# =====================================================

r2cutvec <- 0.05
chrS <- 1     # starting chromosome index (1–22)
plink_path <- "plink"

# Populations to loop through
popvec <- c("EUR", "AFR", "EAS")

refAllelefile <- file.path(out_dir, "reference_alleles.txt")

for (pop in popvec) {

  # set correct reference genotype per population
  if (pop == "EUR") dataname <- bfile_EUR
  if (pop == "AFR") dataname <- bfile_AFR
  if (pop == "EAS") dataname <- bfile_EAS

  cat("\n-----------------------------\n")
  cat("Running LD calculation for:", pop, "\n")
  cat("-----------------------------\n")

  for (rr in 1:length(r2cutvec)) {
    r2cut <- r2cutvec[rr]

    for (chr in chrS:22) {
      # Output name prefix
      outname <- file.path(out_dir, paste0("clumped_chr", chr, "_r2cut", r2cut, "_", pop))

      # Input SNP list for that chromosome and population
      if (r2cut == 1) {
        usedsnpfile <- file.path(out_dir, paste0("extractSNP_chr", chr))
      } else {
        usedsnpfile <- file.path(out_dir, paste0("Chr", chr, "r2cut", r2cut, pop, ".txt"))
      }

      # skip if SNP list file doesn't exist
      if (!file.exists(usedsnpfile)) {
        cat("⚠️ Missing SNP file:", usedsnpfile, "\n")
        next
      }

      # Build PLINK command
      cmd_ld <- paste0(
        plink_path,
        " --bfile ", dataname,
        " --extract ", usedsnpfile,
        " --a1-allele ", refAllelefile,
        " --r2",
        " --ld-window-kb 10000",
        " --ld-window 99999",
        " --ld-window-r2 0",
        " --out ", outname
      )

      cat("Running:", cmd_ld, "\n")
      system(cmd_ld)
    }
  }
}





