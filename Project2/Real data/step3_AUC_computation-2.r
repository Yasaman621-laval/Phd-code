#!/usr/bin/env Rscript

## =================== CLI / Setup ===================
suppressPackageStartupMessages({
  library(data.table)
  library(mvtnorm)
  library(parallel)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript step3_AUC_computation.R <jset>")
jset <- as.numeric(args[1])

file.rjobs   <- "/lustre06/project/6005709/yatah3/real-data/runR/step3/"
inputs_path  <- file.path(file.rjobs, "input.txt")
inputs       <- read.table(inputs_path, header = TRUE, as.is = TRUE)
runIndex     <- inputs[jset, 1]
ttIndex      <- inputs[jset, 2]

cat(">>> Slurm array job index =", jset,
    "| runIndex =", runIndex,
    "| ttIndex =", ttIndex, "\n")

## =================== Paths & constants ===================
dirOutput      <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/"
densityU_dir   <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/densityU_realdata/"
out_dir_test   <- "/lustre06/project/6005709/yatah3/real-data/out-test/"
refAllele_path <- "/lustre06/project/6005709/yatah3/real-data/out/reference_alleles.txt"

taufile_AFR_EUR <- "/lustre06/project/6005709/yatah3/real-data/AFR_EURTau_info_v2.RData"
taufile_EAS_EUR <- "/lustre06/project/6005709/yatah3/real-data/EAS_EURTau_info_v2.RData"

KG_prefix_EUR <- "/lustre06/project/6005709/yatah3/real-data/EUR1000genomes/1kg_hg38_EUR_allchr"
KG_prefix_EAS <- "/lustre06/project/6005709/yatah3/real-data/EAS1000genomes/1kg_hg38_EAS_allchr"
KG_prefix_AFR <- "/lustre06/project/6005709/yatah3/real-data/AFR1000genomes/1kg_hg38_AFR_allchr"

test_gwas_merged <- file.path(out_dir_test, "final_clumped_crosspop.tsv")

SUMMARYAUC_DIR <- "/lustre06/project/6005709/yatah3/real-data/SummaryAUC"
SO_FILE        <- file.path(SUMMARYAUC_DIR, "getAdjCorrelation.so")

savename    <- "step3_AUC_computation"
output_root <- file.path(dirOutput, savename)
dir.create(output_root, showWarnings = FALSE, recursive = TRUE)
R2output <- file.path(output_root, "R2output")
dir.create(R2output, showWarnings = FALSE, recursive = TRUE)

dirtemp <- "/home/yatah3/scratch/"

## =================== Core settings ===================
popvec        <- c("AFR", "EAS", "EUR")
all_usedtrait <- c("1,3", "2,3")
Zscale <- 1; warmStart <- 1; K <- 2; gindex <- 1; plinkver <- 2
cindex <- 1; Before <- 0
rho_vec <- c(seq(0, 0.9, 0.1), 0.95)
singleStart <- 1; PV <- 2
penalty <- "RealmixLOG"

toupper_nt <- function(x) { x <- as.character(x); toupper(x) }

## =================== Which usedtrait & pop? ===================
usedtrait      <- all_usedtrait[runIndex]   # "1,3" or "2,3"
usedtraitsvec  <- as.numeric(strsplit(usedtrait, ",")[[1]])
mainindex      <- usedtraitsvec[1]          # 1=AFR, 2=EAS
secindex       <- usedtraitsvec[2]          # 3=EUR
popuseY_index  <- mainindex                 # evaluate on this pop
POP            <- popvec[popuseY_index]

## =================== Tau load & pick ===================
if (identical(usedtrait, "1,3")) load(taufile_AFR_EUR) else load(taufile_EAS_EUR)
if (!exists("AbsTauvec")) stop("Tau file missing AbsTauvec.")
tau_vec <- AbsTauvec[ttIndex]
tau_vec <- tau_vec[!is.na(tau_vec)]
if (length(tau_vec) == 0) stop("No valid tau values for ttIndex.")
tauuse <- tau_vec[1]

## =================== Test GWAS & N ===================
gwas_all <- fread(test_gwas_merged, data.table = FALSE)

## Build SNPID if missing: CHR:POS:Other:Effect — must match B1/B2 & refAllele
if (!("SNPID" %in% names(gwas_all))) {
  need_cols <- c("CHR","POS","Other_Allele","Effect_Allele")
  miss <- setdiff(need_cols, names(gwas_all))
  if (length(miss)) stop("GWAS missing columns to build SNPID: ", paste(miss, collapse=", "))
  gwas_all$SNPID <- paste0(gwas_all$CHR, ":", gwas_all$POS, ":", 
                           gwas_all$Other_Allele, ":", gwas_all$Effect_Allele)
}

get_test_case_ctrl_N <- function(df_all, pop_index) {
  ncase <- paste0("Ncase_", pop_index)
  nctrl <- paste0("Nctrl_", pop_index)
  ntot  <- paste0("N_",     pop_index)
  if (all(c(ncase, nctrl) %in% names(df_all))) {
    N1 <- round(mean(df_all[[ncase]], na.rm = TRUE))
    N0 <- round(mean(df_all[[nctrl]], na.rm = TRUE))
  } else if (ntot %in% names(df_all)) {
    Ntot <- round(mean(df_all[[ntot]], na.rm = TRUE))
    N1 <- floor(Ntot / 2); N0 <- Ntot - N1
    message("⚠️ Using 50/50 split from test N_: N1=", N1, " N0=", N0)
  } else {
    stop("Test GWAS missing Ncase/Nctrl (or N_) columns.")
  }
  list(N0 = N0, N1 = N1)
}
testNs <- get_test_case_ctrl_N(gwas_all, popuseY_index)
stopifnot(testNs$N0 > 0, testNs$N1 > 0)

## =================== GWAS DF & precompute Z, tau ===================
get_gwas_df <- function(pop_index, gwas_all) {
  data.frame(
    SNP  = gwas_all$SNPID,                              # coord+allele IDs
    A1   = toupper_nt(gwas_all$Effect_Allele),
    BETA = as.numeric(gwas_all[[paste0("beta_", pop_index)]]),
    P    = as.numeric(gwas_all[[paste0("pval_", pop_index)]]),
    MAF  = as.numeric(gwas_all[[paste0("EAF_",  pop_index)]]),
    INFO = 1,
    stringsAsFactors = FALSE
  )
}
gwas_df <- get_gwas_df(popuseY_index, gwas_all)
gwas_df <- gwas_df[complete.cases(gwas_df), ]

Z_vec <- sign(gwas_df$BETA) * qnorm(1 - gwas_df$P/2)
MAFv  <- pmin(gwas_df$MAF, 1 - gwas_df$MAF)
TAUv  <- sqrt(2 * MAFv * (1 - MAFv) * gwas_df$INFO)

## =================== LD reference preload ===================
if (popuseY_index == 1) KG_prefix <- KG_prefix_AFR
if (popuseY_index == 2) KG_prefix <- KG_prefix_EAS
if (popuseY_index == 3) KG_prefix <- KG_prefix_EUR

kg_bim_path <- paste0(KG_prefix, ".bim")
kg_fam_path <- paste0(KG_prefix, ".fam")
kg_bed_path <- paste0(KG_prefix, ".bed")

if (!file.exists(kg_bim_path) || !file.exists(kg_fam_path) || !file.exists(kg_bed_path)) {
  stop("LD reference missing: ", KG_prefix, " (.bim/.fam/.bed)")
}
kg_bim  <- fread(kg_bim_path, header = FALSE, data.table = FALSE)
kg_famN <- nrow(fread(kg_fam_path, header = FALSE))
colnames(kg_bim)[1:6] <- c("CHR","ID","CM","BP","A1","A2")

## =================== Load .so ONCE ===================
if (!file.exists(SO_FILE)) {
  cFile <- sub("\\.so$", ".c", SO_FILE)
  if (file.exists(cFile)) system(paste("R CMD SHLIB", shQuote(cFile)))
}
if (!file.exists(SO_FILE)) stop("Shared library not found: ", SO_FILE)
if (!is.loaded("getAdjCorrelation")) dyn.load(SO_FILE)

## =================== Fast adj. corr. (cached LD) ===================
get_correlation_adj_cached <- function(prs_df, tau_sub, kg_bim, kg_famN, kg_bed_file, pos_thr = 5e8) {
  if (nrow(prs_df) == 0) return(0)

  idx2 <- match(prs_df$SNP, kg_bim$ID)
  keep <- !is.na(idx2)
  if (!any(keep)) return(0)
  prs_df  <- prs_df[keep, , drop=FALSE]
  tau_sub <- tau_sub[keep]
  idx2    <- idx2[keep]

  flip.idx <- which(prs_df$A1 != kg_bim$A1[idx2])
  if (length(flip.idx)) idx2[flip.idx] <- -idx2[flip.idx]

  NSNP <- nrow(prs_df)
  if (NSNP == 0) return(0)

  out <- .C("getAdjCorrelation",
            plinkBed = as.character(kg_bed_file),
            NumSample = as.integer(kg_famN),
            NumSNP    = as.integer(NSNP),
            idx       = as.integer(idx2),
            chrs      = as.integer(kg_bim$CHR[abs(idx2)]),
            pos       = as.integer(kg_bim$BP [abs(idx2)]),
            pos_thr   = as.integer(pos_thr),
            beta_tau  = as.double(prs_df$BETA * tau_sub),
            adj_cor   = as.double(0.1))
  out$adj_cor
}

## =================== AUC (summary) helper ===================
auc_summary_df <- function(prs_df, gwas_df, Z_vec, TAUv, N0, N1,
                           kg_bim, kg_famN, kg_bed_file, pos_thr = 5e8) {
  m <- match(prs_df$SNP, gwas_df$SNP)
  keep <- !is.na(m)
  if (!any(keep)) return(c(NA_real_, NA_real_))
  prs_sub <- prs_df[keep, , drop=FALSE]
  idx_g   <- m[keep]

  gA1  <- gwas_df$A1[idx_g]
  beta <- gwas_df$BETA[idx_g]
  flip <- gA1 != prs_sub$A1
  if (length(flip)) beta[flip] <- -beta[flip]

  Z_sub   <- Z_vec[idx_g]
  tau_sub <- TAUv[idx_g]

  adj_cor <- get_correlation_adj_cached(prs_sub, tau_sub, kg_bim, kg_famN,
                                        kg_bed_file = kg_bed_file, pos_thr = pos_thr)

  b   <- prs_sub$BETA
  num <- sqrt(1/N1 + 1/N0) * sum(b * Z_sub * tau_sub)
  den <- sqrt(2 * sum((b * tau_sub)^2) + 4 * adj_cor)
  if (!is.finite(den) || den == 0) return(c(NA_real_, NA_real_))

  delta <- num / den
  auc0  <- pnorm(delta)
  pxy   <- pmvnorm(lower = c(-delta, -delta), upper = Inf,
                   mean = c(0, 0),
                   sigma = matrix(c(1, 0.5, 0.5, 1), 2, 2))[1]
  var0  <- (pxy - auc0^2) * (N1 + N0) / (N1 * N0)
  c(auc0, var0)
}

## =================== Read refAllele ONCE ===================
refAllele <- fread(refAllele_path, header = FALSE, data.table = FALSE)
colnames(refAllele) <- c("Index","SNP","A1")
# fast lookup table
ref_map <- refAllele[, c("SNP","A1")]
rownames(ref_map) <- ref_map$SNP

## =================== Driver (optimized) ===================
GenPreAUC_Chrs_saveScore_fast <- function(AllBetaMatrix_list,
                                          savef2_list,
                                          ref_map,
                                          gwas_df, Z_vec, TAUv,
                                          N0, N1,
                                          kg_bim, kg_famN, kg_bed_file,
                                          mc.cores = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))) {

  stopifnot(length(AllBetaMatrix_list) == length(savef2_list))

  for (ll in seq_along(AllBetaMatrix_list)) {
    OriBeta <- AllBetaMatrix_list[[ll]]          # rows=SNPIDs (rownames), cols=models (k)
    stopifnot(nrow(OriBeta) > 0, ncol(OriBeta) > 0, !is.null(rownames(OriBeta)))
    AllSNPs <- rownames(OriBeta)
    K       <- ncol(OriBeta)

    savef2 <- savef2_list[[ll]]
    if (file.exists(savef2)) {
      load(savef2)
      if (!exists("AUCsummary") || ncol(AUCsummary) != K) {
        AUCsummary  <- matrix(NA_real_, nrow = 2, ncol = K,
                              dimnames = list(c("AUC", "Var"), paste0("k", seq_len(K))))
        numbetasvec <- rep(NA_integer_, K)
      }
    } else {
      AUCsummary  <- matrix(NA_real_, nrow = 2, ncol = K,
                            dimnames = list(c("AUC", "Var"), paste0("k", seq_len(K))))
      numbetasvec <- rep(NA_integer_, K)
    }

    work_one_k <- function(k) {
      if (!is.na(AUCsummary[1, k])) return(list(k=k, auc=AUCsummary[1,k], var=AUCsummary[2,k], nb=numbetasvec[k]))

      betas_k <- OriBeta[, k]
      if (all(betas_k == 0) || anyNA(betas_k)) return(list(k=k, auc=NA_real_, var=NA_real_, nb=0L))

      snps_k <- AllSNPs
      hasRef <- snps_k %in% ref_map$SNP
      if (!any(hasRef)) return(list(k=k, auc=NA_real_, var=NA_real_, nb=0L))

      snps_k  <- snps_k[hasRef]
      betas_k <- betas_k[hasRef]
      keep    <- betas_k != 0
      if (!any(keep)) return(list(k=k, auc=NA_real_, var=NA_real_, nb=0L))

      a1_ref <- toupper_nt(ref_map[snps_k, "A1", drop = TRUE])

      prs_df <- data.frame(
        SNP  = snps_k[keep],
        A1   = a1_ref[keep],
        BETA = as.numeric(betas_k[keep]),
        stringsAsFactors = FALSE
      )

      res <- tryCatch(
        auc_summary_df(prs_df, gwas_df, Z_vec, TAUv, N0, N1, kg_bim, kg_famN, kg_bed_file),
        error = function(e) c(NA_real_, NA_real_)
      )

      list(k=k, auc=as.numeric(res[1]), var=as.numeric(res[2]), nb=nrow(prs_df))
    }

    res_list <- mclapply(seq_len(K), work_one_k, mc.cores = mc.cores)
    for (r in res_list) {
      AUCsummary[1, r$k] <- r$auc
      AUCsummary[2, r$k] <- r$var
      numbetasvec[r$k]   <- r$nb
    }

    save(AUCsummary, numbetasvec, file = savef2)
    message("✅ Saved results for WB", ll, " → ", savef2)

    rm(OriBeta, AUCsummary, numbetasvec, res_list); gc()
  }
}

## =================== Main loop over rho ===================
ut <- usedtrait
message("==== usedtrait = ", ut, " (", POP, ") ====")

for (rr in seq_along(rho_vec)) {
  rho <- rho_vec[rr]
  message(sprintf(">> tau=%s | rho=%s", as.character(tauuse), as.character(rho)))

  savef2_1 <- file.path(R2output, paste0(
    penalty, "usedtrait", ut, "warmStart", warmStart, "Zscale", Zscale,
    "singleStart", singleStart, "tauuse", tauuse, "rho", rho, "_wB1AUC.RData"
  ))
  savef2_2 <- file.path(R2output, paste0(
    penalty, "usedtrait", ut, "warmStart", warmStart, "Zscale", Zscale,
    "singleStart", singleStart, "tauuse", tauuse, "rho", rho, "_wB2AUC.RData"
  ))
  savef2_list <- list(savef2_1, savef2_2)

  AllBetaMatrix1 <- NULL
  AllBetaMatrix2 <- NULL

  for (chr in 1:22) {
    subtau_saveoutfile <- file.path(
      densityU_dir,
      paste0(penalty, "chr", chr, "usedtrait", ut, "warmStart", warmStart,
             "Zscale", Zscale, "singleStart", singleStart,
             "tauuse", tauuse, "_DensityU.RData")
    )
    if (!file.exists(subtau_saveoutfile)) {
      message("Missing per-chr file: ", subtau_saveoutfile)
      next
    }
    load(subtau_saveoutfile) # provides Allrho_WBMatrix_list
    if (!is.list(Allrho_WBMatrix_list) || length(Allrho_WBMatrix_list) < rr)
      stop("Allrho_WBMatrix_list structure not as expected in: ", subtau_saveoutfile)

    B1 <- Allrho_WBMatrix_list[[rr]][[1]]
    B2 <- Allrho_WBMatrix_list[[rr]][[2]]

    if (is.null(rownames(B1)) || is.null(rownames(B2)))
      stop("No SNP IDs (rownames) in ", subtau_saveoutfile)

    if (!is.null(AllBetaMatrix1) && ncol(B1) != ncol(AllBetaMatrix1))
      stop("Model count differs across chromosomes for WB1 at rho index ", rr, ".")
    if (!is.null(AllBetaMatrix2) && ncol(B2) != ncol(AllBetaMatrix2))
      stop("Model count differs across chromosomes for WB2 at rho index ", rr, ".")

    if (is.null(AllBetaMatrix1)) {
      AllBetaMatrix1 <- B1
      AllBetaMatrix2 <- B2
    } else {
      AllBetaMatrix1 <- rbind(AllBetaMatrix1, B1)
      AllBetaMatrix2 <- rbind(AllBetaMatrix2, B2)
    }

    rm(Allrho_WBMatrix_list, B1, B2)
  }

  if (is.null(AllBetaMatrix1) || is.null(AllBetaMatrix2)) {
    message("No beta matrices assembled; skipping tau/rho combo.")
    next
  }

  GenPreAUC_Chrs_saveScore_fast(
    AllBetaMatrix_list = list(AllBetaMatrix1, AllBetaMatrix2),
    savef2_list        = savef2_list,
    ref_map            = ref_map,
    gwas_df            = gwas_df,
    Z_vec              = Z_vec,
    TAUv               = TAUv,
    N0                 = testNs$N0,
    N1                 = testNs$N1,
    kg_bim             = kg_bim,
    kg_famN            = kg_famN,
    kg_bed_file        = kg_bed_path,
    mc.cores           = max(1L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
  )

  rm(AllBetaMatrix1, AllBetaMatrix2); gc()
}

message("✅ Done: AUC results in ", R2output)
