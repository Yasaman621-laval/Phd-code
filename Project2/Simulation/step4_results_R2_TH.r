
#!/usr/bin/env Rscript
rm(list = ls())

############################################################
# STEP 4 – Collect WB1 & WB2 Results and Plot
# Full structure aligned with PRS-CSx step4 style
# Outputs saved inside: <scenario_dir>/output_step4_WB/
############################################################

########## Libraries ##########
library(dplyr)
library(ggplot2)

########## Read SLURM job input ##########
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing SLURM array index (jset)")

jset <- as.numeric(args[1])

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project2/step4/"
inputs     <- read.table(
  paste0(file.rjobs, "input.txt"),
  header = TRUE,
  as.is  = TRUE
)

if (ncol(inputs) < 2) {
  stop("input.txt must have at least 2 columns: runIndex, scenarioIndex")
}

runIndex      <- inputs[jset, 1]
scenarioIndex <- inputs[jset, 2]

cat(">>> runIndex =", runIndex,
    "| scenarioIndex =", scenarioIndex, "\n")

########## Scenario Folder (same logic as PRS-CSx) ##########
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/"

scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

if (length(scenario_dirs) == 0) {
  stop("No scenario directories found under: ", dirBase)
}
if (scenarioIndex > length(scenario_dirs)) {
  stop("scenarioIndex (", scenarioIndex,
       ") exceeds number of available scenarios (", length(scenario_dirs), ")")
}

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)

cat(">>> Using scenario:", sc_name, "\n")

########## Output Folder ##########
output_step4 <- file.path(sc_dir, "output_step4_WB")
system(paste("mkdir -p", output_step4))

cat(">>> output_step4 WB dir:", output_step4, "\n")

########## R2 output folder (WB1/WB2) ##########
# IMPORTANT: avoid leading '/' inside file.path
R2output <- file.path(sc_dir, "R2output", "R2output1")

cat(">>> R2output directory:", R2output, "\n")

########## General Setup ##########
Zscale   <- 1
savename <- "Trans2"
penalty  <- "RealmixLOG"

popvec          <- c("AFR", "EAS", "EUR")
ordersequse_vec <- c(1, 1)
all_usedtrait   <- c("1,3", "2,3")   # AFR–EUR or EAS–EUR

rho_vec     <- c(seq(0, 0.9, 0.1), 0.95)
singleStart <- 1
warmStart   <- 1

########## Tau File ##########
usedtrait      <- all_usedtrait[runIndex]
usedtraitsvec  <- unlist(strsplit(usedtrait, ","))
usedtraitIndex <- as.numeric(usedtraitsvec)

mainindex <- usedtraitIndex[1]  # 1=AFR or 2=EAS
iiIndexT  <- usedtraitIndex[2]  # 3=EUR

Taufile <- file.path(
  sc_dir,
  paste0(popvec[mainindex], "_", popvec[iiIndexT], "Tau_info_v2.RData")
)

if (!file.exists(Taufile)) {
  stop("Taufile missing: ", Taufile)
}

cat(">>> Loading Taufile:", Taufile, "\n")
load(Taufile)   # must load AbsTauvec (matrix: sim × tt)

if (!exists("AbsTauvec")) {
  stop("AbsTauvec not found in Taufile: ", Taufile)
}

############################################################
# Utility function to collect WB1 or WB2 results
#   WBtype = 1 ? *_wB1R2.RData
#   WBtype = 2 ? *_wB2R2.RData
############################################################

collect_results <- function(WBtype = 1) {

  cat(">>> Collecting WB", WBtype, "results for runIndex =", runIndex, "\n")

  Start     <- 1
  Allresult <- NULL

  # We assume AbsTauvec has dimensions [sim, tt] = [1..10, 1..5]
  nsim <- nrow(AbsTauvec)
  ntt  <- ncol(AbsTauvec)

  for (sim in 1:nsim) {

    bysim_PreR2       <- NULL
    bysim_rhotau      <- NULL
    bysim_numbetasvec <- NULL
    start_bysim       <- 1

    for (tt in 1:ntt) {

      tauuse <- AbsTauvec[sim, tt]

      for (rho in rho_vec) {

        savef2 <- paste0(
          R2output, "/",
          penalty,
          "usedtrait", usedtrait,
          "warmStart",  warmStart,
          "sim",        sim,
          "Zscale",     Zscale,
          "singleStart", singleStart,
          "tauuse",     tauuse,
          "rho",        rho,
          ifelse(WBtype == 1, "_wB1R2.RData", "_wB2R2.RData")
        )

       if (!file.exists(savef2) || file.info(savef2)$size == 0) {
  next
}

        load(savef2)  # must load PreR2, numbetasvec

        if (!exists("PreR2") || !exists("numbetasvec")) {
          stop("File ", savef2, " does not contain PreR2 and numbetasvec")
        }

        if (length(numbetasvec) != ncol(PreR2)) {
          stop("Length mismatch in ", savef2,
               ": length(numbetasvec)=", length(numbetasvec),
               " but ncol(PreR2)=", ncol(PreR2))
        }


# ---- PRINT R2 FOR DEBUGGING (for each sim, tau, rho) ----

whmax_local <- which.max(PreR2[1, ])   # best validation R²

cat("\n----------------------------\n")
cat("sim =", sim,
    "| tau =", tauuse,
    "| rho =", rho, "\n")
cat("File:", savef2, "\n")

cat("whmax (best column index):", whmax_local, "\n")
cat("R2_row1_best =", PreR2[1, whmax_local], "\n")
cat("R2_row2_best =", PreR2[2, whmax_local], "\n\n")
cat("----------------------------\n\n")


        temprt <- rbind(
          rep(rho,    ncol(PreR2)),  # row1 = rho
          rep(tauuse, ncol(PreR2))   # row2 = tau
        )

        if (start_bysim == 1) {
          bysim_PreR2       <- PreR2
          bysim_rhotau      <- temprt
          bysim_numbetasvec <- numbetasvec
          start_bysim       <- 0
        } else {
          bysim_PreR2       <- cbind(bysim_PreR2,       PreR2)
          bysim_rhotau      <- cbind(bysim_rhotau,      temprt)
          bysim_numbetasvec <- c(bysim_numbetasvec,     numbetasvec)
        }

        rm(PreR2, numbetasvec)
      }
    }

    # If no files found for this sim, skip it
    if (is.null(bysim_PreR2) || length(bysim_PreR2) == 0) {
      cat("   [WB", WBtype, "] No R2 files found for sim =", sim, "\n")
      next
    }

    # bysim_PreR2 has rows: [1] validation R2?, [2] test R2? (keep as original)
    # The original code used row 2 as the final R2; we preserve that.
    whmax <- which.max(bysim_PreR2[1, ])

    temp <- c(
      tau_best     = bysim_rhotau[2, whmax],   # best tau
      rho_best     = bysim_rhotau[2, whmax],   # best rho
      R2_best      = bysim_PreR2[2, whmax],    # R2 (second row)
      nBetas_best  = bysim_numbetasvec[whmax], # num betas at best tuning
      nBetas_max   = max(bysim_numbetasvec),   # max betas across grid
      pop_code     = usedtraitIndex[1]         # 1=AFR or 2=EAS
    )

    if (Start == 1) {
      Allresult <- temp
      Start     <- 0
    } else {
      Allresult <- rbind(Allresult, temp)
    }
  }

  if (is.null(Allresult)) {
    stop("No WB", WBtype, " results found for runIndex=", runIndex,
         " in scenario=", sc_name)
  }

  Allresult <- as.matrix(Allresult)
  colnames(Allresult) <- c(
    "tau_best", "rho_best", "R2_best",
    "nBetas_best", "nBetas_max", "pop_code"
  )

  return(Allresult)
}

############################################################
# Run WB1 and WB2 data collection
############################################################

AllresultB1 <- collect_results(WBtype = 1)
AllresultB2 <- collect_results(WBtype = 2)

save(
  AllresultB1,
  file = file.path(output_step4,
                   paste0("AllresultWB1_runIndex", runIndex, ".RData"))
)
save(
  AllresultB2,
  file = file.path(output_step4,
                   paste0("AllresultWB2_runIndex", runIndex, ".RData"))
)

cat(">>> Summary WB1 (tau, rho, R2, nBetas_best, nBetas_max):\n")
print(apply(AllresultB1[, 1:5, drop = FALSE], 2, summary))

cat(">>> Summary WB2 (tau, rho, R2, nBetas_best, nBetas_max):\n")
print(apply(AllresultB2[, 1:5, drop = FALSE], 2, summary))

############################################################
# Optional: WB1 vs WB2 boxplot of R² (aligned with PRS-CSx style)
############################################################

R2_WB1 <- as.numeric(AllresultB1[, "R2_best"])
R2_WB2 <- as.numeric(AllresultB2[, "R2_best"])

Method <- c(
  rep("WB1", length(R2_WB1)),
  rep("WB2", length(R2_WB2))
)

R2 <- c(R2_WB1, R2_WB2)

df_plot <- data.frame(Method, R2)
df_plot$Method <- factor(df_plot$Method, levels = c("WB1", "WB2"))

p <- ggplot(df_plot, aes(x = Method, y = R2)) +
  geom_boxplot(color = "black", fill = "white",
               outlier.shape = NA, size = 0.8) +
  labs(
    title = paste0("WB R\u00B2 (Scenario: ", sc_name,
                   ", runIndex=", runIndex, ")"),
    x = "Wavelet Block (WB1 vs WB2)",
    y = expression(R^2)
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.line  = element_line(size = 0.8),
    axis.ticks = element_line(size = 0.8),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

out_png <- file.path(output_step4,
                     paste0("WB_R2_boxplot_runIndex", runIndex, ".png"))
ggsave(out_png, p, width = 8, height = 5, dpi = 300)

cat(">>> WB1 vs WB2 R2 boxplot saved to:", out_png, "\n")
cat(">>> STEP4 WB (WB1 + WB2) Completed for scenario =", sc_name,
    "runIndex =", runIndex, "\n")
