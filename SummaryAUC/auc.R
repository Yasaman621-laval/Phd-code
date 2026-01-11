# ================================================================
# Cleaned version of auc.R (Yasaman adjusted)
# Purpose: compute AUC from PRS model & validation GWAS summary
# ================================================================

get.correlation.adj <- function(prs.model,
                                tau,
                                KG.plink.pre,
                                soFile,
                                pos_thr = 5e8) {

  kg.fam <- read.table(file = paste0(KG.plink.pre, ".fam"), header = FALSE, stringsAsFactors = FALSE)
  NSample <- nrow(kg.fam)

  kg.bim <- read.table(file = paste0(KG.plink.pre, ".bim"), header = FALSE, stringsAsFactors = FALSE)
  snps <- intersect(prs.model[, 1], kg.bim[, 2])
  NSNP <- length(snps)
  idx <- match(snps, prs.model[, 1])
  prs.model <- prs.model[idx, ]
  tau <- tau[idx]

  idx <- match(prs.model[, 1], kg.bim[, 2])
  chrs <- kg.bim[idx, 1]
  pos <- kg.bim[idx, 4]
  flip.idx <- which(prs.model[, 2] != kg.bim[idx, 5])
  idx[flip.idx] <- -idx[flip.idx]

  beta_tau <- prs.model[, 3] * tau
  kg.bed.file <- paste0(KG.plink.pre, ".bed")

  # Compile C file if needed
  if (!file.exists(soFile)) {
    cFile <- gsub(".so$", ".c", soFile)
    system(paste("R CMD SHLIB", cFile))
  }
  dyn.load(soFile)

  adj.cor.results <- .C("getAdjCorrelation",
                        plinkBed = kg.bed.file,
                        NumSample = as.integer(NSample),
                        NumSNP = as.integer(NSNP),
                        idx = as.integer(idx),
                        chrs = as.integer(chrs),
                        pos = as.integer(pos),
                        pos_thr = as.integer(pos_thr),
                        beta_tau = as.double(beta_tau),
                        adj_cor = as.double(0.1))

  return(adj.cor.results$adj_cor)
}

# ----------------------------------------------------------------
# Main function: compute AUC and its variance
# ----------------------------------------------------------------
auc <- function(prs.model.file,
                gwas.summary.stats.file,
                N0,
                N1,
                soFile = "getAdjCorrelation.so",
                KG.plink.pre = NULL,
                pos_thr = 5e8,
                flag.correlation.adj.imputated.data = FALSE) {

  library("mvtnorm")

  prs.model <- read.delim(file = prs.model.file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  gwas.summary.stats <- read.table(file = gwas.summary.stats.file, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

  # Align SNPs between PRS model and GWAS summary
  snps <- intersect(prs.model[, 1], gwas.summary.stats[, "SNP"])
  idx <- match(snps, gwas.summary.stats[, "SNP"])
  gwas.summary.stats <- gwas.summary.stats[idx, ]
  idx <- match(snps, prs.model[, 1])
  prs.model <- prs.model[idx, ]

  # Determine beta values
  if ("OR" %in% colnames(gwas.summary.stats)) {
    beta_gwas <- log(gwas.summary.stats[, "OR"])
  } else {
    beta_gwas <- gwas.summary.stats[, "BETA"]
  }

  # Align alleles
  idx <- gwas.summary.stats[, "A1"] != prs.model[, 2]
  beta_gwas[idx] <- -beta_gwas[idx]

  P <- gwas.summary.stats[, "P"]
  MAF <- gwas.summary.stats[, "MAF"]
  Z <- sign(beta_gwas) * qnorm(1 - P / 2)

  if ("INFO" %in% colnames(gwas.summary.stats)) {
    INFO <- as.numeric(gwas.summary.stats[, "INFO"])
  } else {
    INFO <- rep(1, length(MAF))
  }

  tau <- sqrt(2 * MAF * (1 - MAF) * INFO)
  beta_prs <- prs.model[, 3]

  # ---- correlation adjustment ----
  if (flag.correlation.adj.imputated.data) {
    stop("Imputed dosage-based correlation mode not supported in this simplified version. Use LD mode instead.")
  } else {
    if (!is.null(KG.plink.pre)) {
      adj.cor <- get.correlation.adj(prs.model, tau, KG.plink.pre, soFile, pos_thr)
    } else {
      adj.cor <- 0  # independent-SNP assumption
    }
  }

  delta <- sqrt(1 / N1 + 1 / N0) * sum(beta_prs * Z * tau) /
    sqrt(2 * sum((beta_prs * tau)^2) + 4 * adj.cor)

  auc0 <- pnorm(delta)
  pxy <- pmvnorm(lower = c(-delta, -delta),
                 upper = Inf,
                 mean = c(0, 0),
                 sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2))[1]
  var0 <- (pxy - auc0^2) * (N1 + N0) / (N1 * N0)

  return(c(auc0, var0))
}

# ================================================================
# End of file (no demo run below)
# ================================================================
