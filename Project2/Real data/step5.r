rm(list=ls())
# ================================================================
#   Step 1. Load packages
# ================================================================
library(stringr)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)

# ================================================================
#   Step 2. Adjust these settings
# ================================================================
# Path to your R2output folder (same for both AFR–EUR and EAS–EUR)
r2_dir <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/step3_AUC_computation/R2output"

# Change this line depending on which pair you analyze:
usedtrait <- "2,3"   # "1,3" = AFR–EUR   or   "2,3" = EAS–EUR

# ================================================================
#   Step 3. List all RData files for that usedtrait
# ================================================================
pattern <- paste0("usedtrait", usedtrait, ".*_wB[12]AUC\\.RData$")
files <- list.files(r2_dir, pattern = pattern, full.names = TRUE)
if (length(files) == 0) stop("No AUC files found for usedtrait ", usedtrait)

# Extract tau, rho, WB index (1 = within-pop, 2 = EUR→target)
info <- tibble(
  file = files,
  WB   = str_extract(files, "_wB[12]AUC") |> str_extract("[12]"),
  tau  = str_extract(files, "(?<=tauuse)\\d+\\.?\\d*(e-?\\d+)?"),
  rho  = str_extract(files, "(?<=rho)\\d*\\.?\\d*")
) |> mutate(
  tau = as.numeric(tau),
  rho = as.numeric(rho),
  WB  = ifelse(WB == "1", "WB1_within", "WB2_cross")
)

# ================================================================
#   Step 4. Function to read one AUC file
# ================================================================
read_auc <- function(path) {
  env <- new.env()
  load(path, envir = env)
  auc <- if (exists("AUCsummary", envir = env)) env$AUCsummary["AUC", ] else NA_real_
  nb  <- if (exists("numbetasvec", envir = env)) env$numbetasvec else NA_integer_
  tibble(k = seq_along(nb), AUC = auc, NumBetas = nb)
}

# ================================================================
#   Step 5. Combine all results
# ================================================================
all_results <- info %>%
  rowwise() %>%
  mutate(content = list(read_auc(file))) %>%
  unnest(cols = c(content))

summary_tbl <- all_results %>%
  group_by(WB, tau, rho) %>%
  summarise(
    mean_AUC = mean(AUC, na.rm = TRUE),
    median_AUC = median(AUC, na.rm = TRUE),
    total_NumBetas = sum(NumBetas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(WB, tau, rho)

# ================================================================
#   Step 6. Create wide table: WB1 vs WB2
# ================================================================
comparison_tbl <- summary_tbl %>%
  select(WB, tau, rho, total_NumBetas, mean_AUC) %>%
  pivot_wider(
    names_from = WB,
    values_from = c(total_NumBetas, mean_AUC)
  ) %>%
  arrange(tau, rho)

# ================================================================
#   Step 7. Compute differences and ratios
# ================================================================
comparison_tbl <- comparison_tbl %>%
  mutate(
    diff_AUC = mean_AUC_WB2_cross - mean_AUC_WB1_within,
    transfer_ratio = mean_AUC_WB2_cross / mean_AUC_WB1_within,
    diff_NumBetas = total_NumBetas_WB2_cross - total_NumBetas_WB1_within
  )
comparison_tbl <- comparison_tbl %>%
  filter(!is.na(tau) & !is.na(rho))


#==============================================================
#   Step 8. Plot ΔAUC
# ================================================================


# Load libraries
library(ggplot2)

# (Assuming you already have 'comparison_tbl' and 'usedtrait' in your workspace)

# Define output file path (adjust as you like)
output_path <- paste0("AUC_diff_plot_", gsub(",", "_", usedtrait), ".png")

# Create the plot
p <- ggplot(comparison_tbl, aes(x = rho, y = diff_AUC, color = factor(tau))) +
  geom_line(size = 1.1) +
  geom_point(size = 1.5) +
  theme_minimal(base_size = 13) +
  labs(
    title = paste0(
      "ΔAUC (EUR→", ifelse(usedtrait == "1,3", "AFR", "EAS"),
      " minus ", ifelse(usedtrait == "1,3", "AFR→AFR", "EAS→EAS"), ")"
    ),
    x = expression(rho),
    y = expression(Delta * "AUC (cross - within)"),
    color = expression(tau)
  )

# Save as high-resolution PNG
ggsave(
  filename = output_path,
  plot = p,
  width = 8, height = 6, units = "in",
  dpi = 300
)

cat("✅ Plot saved to:", output_path, "\n")

# ================================================================
#   Step 9. Plot transferability ratio (AUC_cross / AUC_within)
# ================================================================
output_path2 <- paste0("transferability ratio_plot_", gsub(",", "_", usedtrait), ".png")


g = ggplot(comparison_tbl, aes(x = rho, y = transfer_ratio, color = factor(tau))) +
  geom_line(size = 1.1) +
  geom_point(size = 1.5) +
  theme_minimal(base_size = 13) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray40") +
  labs(
    title = paste0("Transferability ratio (AUCcross / AUCwithin): ",
                   ifelse(usedtrait == "1,3", "AFR–EUR", "EAS–EUR")),
    x = expression(rho),
    y = "Transferability ratio",
    color = expression(tau)
  )
ggsave(
  filename = output_path2,
  plot = g,
  width = 8, height = 6, units = "in",
  dpi = 300
)

# ================================================================
#   Step 10. Print numeric summaries
# ================================================================
cat("\n===== Mean ΔAUC per τ =====\n")
print(
  comparison_tbl %>%
    group_by(tau) %>%
    summarise(mean_diff = mean(diff_AUC, na.rm = TRUE),
              mean_ratio = mean(transfer_ratio, na.rm = TRUE))
)

cat("\n✅ Finished summary for usedtrait =", usedtrait, "\n")
