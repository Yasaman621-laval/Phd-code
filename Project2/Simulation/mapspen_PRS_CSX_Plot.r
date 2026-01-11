#!/usr/bin/env Rscript
rm(list = ls())

library(dplyr)
library(ggplot2)

############################################################
# INPUT DIRECTORIES
############################################################

scenario_dir <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/sim_hsq0.9_rho0.9_train25000-25000-200000"
output_step4_mapspen <- paste0(scenario_dir, "/output_step4_WB/")
output_step4_prscsx  <- paste0(scenario_dir, "/output_step4_PRScsx/")

############################################################
# 1) LOAD MAPSPEN RESULTS
############################################################

load(paste0(output_step4_mapspen, "AllresultWB1_runIndex1.RData"))   # AFR WB1
MAPSPEN_AFR       <- AllresultB1[,3]

load(paste0(output_step4_mapspen, "AllresultWB2_runIndex1.RData"))   # EUR?AFR WB2
MAPSPEN_EUR_AFR   <- AllresultB2[,3]

load(paste0(output_step4_mapspen, "AllresultWB1_runIndex2.RData"))   # EAS WB1
MAPSPEN_EAS       <- AllresultB1[,3]

load(paste0(output_step4_mapspen, "AllresultWB2_runIndex2.RData"))   # EUR?EAS WB2
MAPSPEN_EUR_EAS   <- AllresultB2[,3]

############################################################
# 2) LOAD PRS-CSx RESULTS
############################################################

load(paste0(output_step4_prscsx, "Allresult_PRScsx_AFR_within_nonmeta.RData"))
load(paste0(output_step4_prscsx, "Allresult_PRScsx_EUR_to_AFR_nonmeta.RData"))
load(paste0(output_step4_prscsx, "Allresult_PRScsx_AFR_meta.RData"))

AFR_nonmeta  <- AFR_within_R2
AFR_cross    <- EUR_to_AFR_R2
AFR_meta     <- AFR_meta_R2

load(paste0(output_step4_prscsx, "Allresult_PRScsx_EAS_within_nonmeta.RData"))
load(paste0(output_step4_prscsx, "Allresult_PRScsx_EUR_to_EAS_nonmeta.RData"))
load(paste0(output_step4_prscsx, "Allresult_PRScsx_EAS_meta.RData"))

EAS_nonmeta  <- EAS_within_R2
EAS_cross    <- EUR_to_EAS_R2
EAS_meta     <- EAS_meta_R2

############################################################
# 3) BUILD DATA FRAME
############################################################

Population <- c(
  # AFR block
  rep("AFR (MAPSPEN)",      length(MAPSPEN_AFR)),
  rep("AFR (PRS-CSx)",      length(AFR_nonmeta)),
  rep("AFR (PRS-CSx meta)", length(AFR_meta)),

  # EUR?AFR block
  rep("EUR-AFR (MAPSPEN)",  length(MAPSPEN_EUR_AFR)),
  rep("EUR-AFR (PRS-CSx)",  length(AFR_cross)),

  # EAS block
  rep("EAS (MAPSPEN)",      length(MAPSPEN_EAS)),
  rep("EAS (PRS-CSx)",      length(EAS_nonmeta)),
  rep("EAS (PRS-CSx meta)", length(EAS_meta)),

  # EUR?EAS block
  rep("EUR-EAS (MAPSPEN)",  length(MAPSPEN_EUR_EAS)),
  rep("EUR-EAS (PRS-CSx)",  length(EAS_cross))
)

R2 <- c(
  MAPSPEN_AFR,
  AFR_nonmeta,
  AFR_meta,

  MAPSPEN_EUR_AFR,
  AFR_cross,

  MAPSPEN_EAS,
  EAS_nonmeta,
  EAS_meta,

  MAPSPEN_EUR_EAS,
  EAS_cross
)

Method <- c(
  # AFR block
  rep("MAPSPEN",      length(MAPSPEN_AFR)),
  rep("PRS-CSx",      length(AFR_nonmeta)),
  rep("PRS-CSx-meta", length(AFR_meta)),

  # EUR?AFR block
  rep("MAPSPEN",      length(MAPSPEN_EUR_AFR)),
  rep("PRS-CSx",      length(AFR_cross)),

  # EAS block
  rep("MAPSPEN",      length(MAPSPEN_EAS)),
  rep("PRS-CSx",      length(EAS_nonmeta)),
  rep("PRS-CSx-meta", length(EAS_meta)),

  # EUR?EAS block
  rep("MAPSPEN",      length(MAPSPEN_EUR_EAS)),
  rep("PRS-CSx",      length(EAS_cross))
)

df <- data.frame(Population, Method, R2)

df$Population <- factor(df$Population, levels = c(
  "AFR (MAPSPEN)",
  "AFR (PRS-CSx)",
  "AFR (PRS-CSx meta)",
  "EUR-AFR (MAPSPEN)",
  "EUR-AFR (PRS-CSx)",
  "EAS (MAPSPEN)",
  "EAS (PRS-CSx)",
  "EAS (PRS-CSx meta)",
  "EUR-EAS (MAPSPEN)",
  "EUR-EAS (PRS-CSx)"
))

############################################################
# 4) MEAN LABELS
############################################################

df_mean <- df %>%
  group_by(Population) %>%
  summarise(mean_R2 = mean(R2))

############################################################
# 5) FINAL PLOT
############################################################

p <- ggplot(df, aes(x = Population, y = R2, fill = Method)) +
  geom_boxplot(outlier.shape = NA, width = 0.65) +
  geom_text(
    data = df_mean,
    aes(x = Population, y = mean_R2, label = sprintf("%.3f", mean_R2)),
    inherit.aes = FALSE,
    vjust = -1.2,
    size = 3.8,
    fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "MAPSPEN"      = "#1f77b4",
    "PRS-CSx"      = "#B8860B",
    "PRS-CSx-meta" = "coral"
  )) +
  labs(
    title = "Comparison of MAPSPEN vs PRS-CSx Across Populations",
    x = "",
    y = expression(R^2),
    fill = "Method"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

############################################################
# 6) SAVE
############################################################

ggsave(
  "/lustre09/project/6005709/yatah3/simulation/project2/step5_finalplot/MAPSPEN_vs_PRScsx_scenario3.png",
  p, width = 14, height = 7, dpi = 300
)

cat(">>> Final MAPSPEN vs PRS-CSx comparison saved.\n")

