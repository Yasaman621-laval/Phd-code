rm(list = ls())

# ============================================================
# MAPSPEN vs PRS-CSx — R² distributions (REAL DATA)
# MAPSPEN: varies over (tau, rho)
# PRS-CSx: varies over phi
# ============================================================

suppressPackageStartupMessages({
  library(stringr)
  library(dplyr)
  library(purrr)
  library(ggplot2)
})

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
dir_mapspen <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/step3_AUC_computation/R2output"
dir_prscsx  <- "/lustre03/project/6005709/yatah3/real-data/step3_AUC_computation/R2output-PRCCSX"

# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
safe_load <- function(f) {
  e <- new.env(parent = emptyenv())
  ok <- tryCatch({ load(f, envir = e); TRUE }, error = function(x) FALSE)
  if (!ok) return(NULL)
  as.list(e)
}

# best R² (AUC row) WITHIN ONE FILE
get_best_R2 <- function(obj) {
  if (!"AUCsummary" %in% names(obj)) return(NA_real_)
  x <- obj$AUCsummary
  if (!("AUC" %in% rownames(x))) return(NA_real_)
  max(as.numeric(x["AUC", ]), na.rm = TRUE)
}

# ------------------------------------------------------------
# MAPSPEN — keep ALL (tau, rho)
# ------------------------------------------------------------
files_mapspen <- list.files(
  dir_mapspen,
  pattern = "RealmixLOGusedtrait[12],3.*_wB[12]AUC\\.RData$",
  full.names = TRUE
)

mapspen_df <- map_dfr(files_mapspen, function(f) {

  obj <- safe_load(f)
  if (is.null(obj)) return(NULL)

  val <- get_best_R2(obj)
  if (!is.finite(val)) return(NULL)

  bn <- basename(f)

  # Extract tau, rho
  tau <- as.numeric(str_extract(bn, "(?<=tauuse)[0-9.eE-]+"))
  rho <- as.numeric(str_extract(bn, "(?<=rho)[0-9.]+"))

  # Trait pairing
  trait_pair <- case_when(
    str_detect(bn, "usedtrait1,3") ~ "AFR_EUR",
    str_detect(bn, "usedtrait2,3") ~ "EAS_EUR",
    TRUE ~ NA_character_
  )

  wb <- ifelse(str_detect(bn, "_wB1"), "WB1", "WB2")

  group <- case_when(
    trait_pair == "AFR_EUR" & wb == "WB1" ~ "AFR",
    trait_pair == "AFR_EUR" & wb == "WB2" ~ "EUR-AFR",
    trait_pair == "EAS_EUR" & wb == "WB1" ~ "EAS",
    trait_pair == "EAS_EUR" & wb == "WB2" ~ "EUR-EAS",
    TRUE ~ NA_character_
  )

  tibble(
    method = "MAPSPEN",
    group  = group,
    tau    = tau,
    rho    = rho,
    value  = val
  )
})

# ------------------------------------------------------------
# PRS-CSx — keep ALL phi
# ------------------------------------------------------------
files_prscsx <- list.files(
  dir_prscsx,
  pattern = "^PRSCSx_.*_AUC\\.RData$",
  full.names = TRUE
)

prscsx_df <- map_dfr(files_prscsx, function(f) {

  obj <- safe_load(f)
  if (is.null(obj)) return(NULL)

  val <- get_best_R2(obj)
  if (!is.finite(val)) return(NULL)

  bn <- basename(f)

  phi <- as.numeric(str_extract(bn, "(?<=phiuse)[0-9.eE-]+"))
  is_meta <- str_detect(bn, "_META_")

  group <- case_when(
    str_detect(bn, "AFR-AFR") ~ "AFR",
    str_detect(bn, "EAS-EAS") ~ "EAS",
    str_detect(bn, "AFR-EUR") ~ "EUR-AFR",
    str_detect(bn, "EAS-EUR") ~ "EUR-EAS",
    str_detect(bn, "AFR_META") ~ "AFR",
    str_detect(bn, "EAS_META") ~ "EAS",
    TRUE ~ NA_character_
  )

  tibble(
    method = ifelse(is_meta, "PRS-CSx-meta", "PRS-CSx"),
    group  = group,
    phi    = phi,
    value  = val
  )
})

# ------------------------------------------------------------
# Combine (NO COLLAPSING)
# ------------------------------------------------------------
all_df <- bind_rows(mapspen_df, prscsx_df) %>%
  filter(!is.na(group))

# ------------------------------------------------------------
# Exact x-axis order (matches your figure)
# ------------------------------------------------------------
x_order <- c(
  "AFR (MAPSPEN)", "AFR (PRS-CSx)", "AFR (PRS-CSx-meta)",
  "EUR-AFR (MAPSPEN)", "EUR-AFR (PRS-CSx)",
  "EAS (MAPSPEN)", "EAS (PRS-CSx)", "EAS (PRS-CSx-meta)",
  "EUR-EAS (MAPSPEN)", "EUR-EAS (PRS-CSx)"
)

all_df <- all_df %>%
  mutate(
    xcat = factor(paste0(group, " (", method, ")"), levels = x_order),
    method = factor(method, levels = c("MAPSPEN", "PRS-CSx", "PRS-CSx-meta"))
  ) %>%
  filter(!is.na(xcat))

# ------------------------------------------------------------
# Median labels
# ------------------------------------------------------------
lab_df <- all_df %>%
  group_by(xcat) %>%
  summarise(
    med = median(value),
    q3  = quantile(value, 0.75),
    .groups = "drop"
  )

# ------------------------------------------------------------
# Group separators
# ------------------------------------------------------------
separators <- data.frame(x = c(3.5, 5.5, 8.5))

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
p <- ggplot(all_df, aes(x = xcat, y = value, fill = method)) +

  geom_boxplot(width = 0.7, outlier.shape = NA) +

  geom_text(
    data = lab_df,
    aes(x = xcat, y = q3, label = sprintf("%.3f", med)),
    inherit.aes = FALSE,
    vjust = -0.4,
    size = 4,
    fontface = "bold"
  ) +

  geom_vline(
    data = separators,
    aes(xintercept = x),
    linewidth = 0.6,
    color = "black"
  ) +

  labs(
    title = "Comparison of MAPSPEN vs PRS-CSx Across Populations",
    x = NULL,
    y = expression(AUC),
    fill = "Method"
  ) +

  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )


 ggsave("MAPSPEN_vs_PRSCSx_AUC_DISTRIBUTION.png", p, width = 16, height = 6, dpi = 300)
