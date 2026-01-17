rm(list=ls())
library(readr)
library(dplyr)
library(ggplot2)

# ===============================
# 1) Load GWAS files
# ===============================
eur <- read_tsv("eur-t2d-train.h.tsv.gz")
eas <- read_tsv("eas-t2d-train.h.tsv.gz")
afr <- read_tsv("afr-t2d-train.h.tsv.gz")

# ===============================
# 2) Harmonize EUR (already reported beta, but might be scale-mismatched)
# ===============================
eur <- eur %>%
  mutate(
    # Step 1: Z from p-value (scale-free)
    Z_from_pval = sign(beta) * qnorm(1 - p_value / 2),

    # Step 2: recalc beta on log(OR) scale
    beta_logOR = Z_from_pval / sqrt(2 * N * effect_allele_frequency * (1 - effect_allele_frequency)),

    # Step 3: recompute SE
    se_recalc   = 1 / sqrt(2 * N * effect_allele_frequency * (1 - effect_allele_frequency)),

    # Step 4: sanity-check p-values
    pval_check = 2 * pnorm(-abs(beta_logOR / se_recalc))
  )

# ===============================
# 3) Harmonize EAS (odds ratios ? log(OR))
# ===============================
eas <- eas %>%
  mutate(
    beta_logOR = log(odds_ratio),
    Z_from_pval = sign(beta_logOR) * qnorm(1 - p_value / 2),
    se_recalc   = abs(beta_logOR / Z_from_pval),  # back-calc SE
    pval_check  = 2 * pnorm(-abs(Z_from_pval))
  )

# ===============================
# 4) Harmonize AFR (odds ratios ? log(OR))
# ===============================
afr <- afr %>%
  mutate(
    beta_logOR = log(odds_ratio),
    Z_from_pval = sign(beta_logOR) * qnorm(1 - p_value / 2),
    se_recalc   = abs(beta_logOR / Z_from_pval),
    pval_check  = 2 * pnorm(-abs(Z_from_pval))
  )

# ===============================
# 5) QC checks
# ===============================
# EUR
print(cor(eur$p_value, eur$pval_check, use="complete.obs"))
summary(eur$beta_logOR)
hist(eur$beta_logOR, breaks=100, main="EUR recalculated betas")

# EAS
print(cor(eas$p_value, eas$pval_check, use="complete.obs"))
summary(eas$beta_logOR)
hist(eas$beta_logOR, breaks=100, main="EAS recalculated betas")

# AFR
print(cor(afr$p_value, afr$pval_check, use="complete.obs"))
summary(afr$beta_logOR)
hist(afr$beta_logOR, breaks=100, main="AFR recalculated betas")

# ===============================
# 6) Save harmonized outputs
# ===============================
write_tsv(eur, "EUR_harmonized.tsv")
write_tsv(eas, "EAS_harmonized.tsv")
write_tsv(afr, "AFR_harmonized.tsv")



#

