rm(list = ls())
library(tidyr)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)

# ------------------------------------------------------------
# 1) Read SLURM job input
# ------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) stop("Missing SLURM array index")
jset <- as.numeric(args[1])

file.rjobs <- "/lustre09/project/6005709/yatah3/simulation/project1/step5-plot/"
inputs     <- read.table(paste0(file.rjobs, "input.txt"), header = TRUE)

scenarioIndex <- inputs[jset, 1]

# ------------------------------------------------------------
# 2) Identify scenario directory
# ------------------------------------------------------------
dirBase <- "/lustre09/project/6005709/yatah3/simulation/SimuGenotype/three_traits/"
scenario_dirs <- list.dirs(dirBase, full.names = TRUE, recursive = FALSE)
scenario_dirs <- scenario_dirs[grepl("sim_hsq", basename(scenario_dirs))]

sc_dir  <- scenario_dirs[scenarioIndex]
sc_name <- basename(sc_dir)

cat(">>> Scenario:", sc_name, "\n")

# ------------------------------------------------------------
# 3) Extract true hsq and rho from folder name
# ------------------------------------------------------------
hsq <- as.numeric(str_extract(sc_name, "(?<=hsq)[0-9.]+"))
rho <- as.numeric(str_extract(sc_name, "(?<=rho)[0-9.]+"))

cat(">>> hsq =", hsq, " rho =", rho, "\n")

# ------------------------------------------------------------
# 4) Load TRUE p-values for this scenario
# ------------------------------------------------------------
df <- read.csv(paste0(dirBase, "SimulationSummary_threeTraits.csv"))

get_pi_values <- function(df, hsq_value, rho_value) {
  row <- df[df$hsq == hsq_value & df$rho == rho_value, ]
  if (nrow(row) == 0) stop("No match found in SimulationSummary.")
  row[, c("pi_0","pi_1","pi_2","pi_3","pi_12","pi_13","pi_23","pi_123")]
}

true_pi_wide <- get_pi_values(df, hsq, rho)

pi_values <- data.frame(
  Parameter = names(true_pi_wide),
  Real_Value = as.numeric(true_pi_wide[1, ])
)

cat(">>> True p-values:\n")
print(pi_values)

# ------------------------------------------------------------
# 5) Load es-alpha results for sim = 1..10
# ------------------------------------------------------------
sim_list <- list()

for (sim in 1:10) {

  pattern <- paste0("es_alpha_best_sim", sim,
                    "_scenario", scenarioIndex, "_row")

  files <- list.files(sc_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    cat(">>> Missing es-alpha for sim", sim, "\n")
    next
  }

  load(files[1])  # loads: es_alphaMatrix, row_index, tauuse

  selected_row <- es_alphaMatrix[row_index, ]

  sim_list[[sim]] <- selected_row
}

# Convert list ? matrix ? long table
sim_mat <- do.call(rbind, sim_list)
colnames(sim_mat) <- c("pi_0","pi_1","pi_2","pi_3","pi_12","pi_13","pi_23","pi_123")

sim_df <- cbind(sim = 1:nrow(sim_mat), sim_mat)

sim_df <- as.data.frame(sim_df)


data_long <- pivot_longer(sim_df,
                          cols = starts_with("pi_"),
                          names_to = "Parameter",
                          values_to = "Value")

# ------------------------------------------------------------
# 6) Mean per parameter
# ------------------------------------------------------------
sim_means <- data_long %>%
  group_by(Parameter) %>%
  summarise(Mean_Value = mean(Value), .groups = "drop")

# ------------------------------------------------------------
# 7) JASA-style boxplot
# ------------------------------------------------------------
p <- ggplot(data_long, aes(Parameter, Value)) +
  geom_boxplot(fill = "white", color = "black",outlier.shape = NA) +

  # Add means with mapped color for legend
  geom_point(
    data = sim_means,
    aes(Parameter, Mean_Value, color = "Estimated Mean"),
    size = 3
  ) +

  # Add true values with mapped color for legend
  geom_point(
    data = pi_values,
    aes(Parameter, Real_Value, color = "True Value"),
    size = 3
  ) +

  # Define legend colors & labels
  scale_color_manual(
    name = "Legend",
    values = c(
      "Estimated Mean" = "darkgoldenrod3",
      "True Value"     = "navy"
    )
  ) +

  theme_classic(base_size = 14) +
  labs(
    x = "p-pattern",
    y = "Posterior Inclusion Probability",
   title = paste0(
  "Estimated p-pattern distribution across simulations (Scenario ", scenarioIndex, ")\n",
  "hsq = ", hsq, ",  rho = ", rho, "  —  Three traits"
)
  ) +
  theme(
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    legend.position = c(0.9,0.8),
     legend.title = element_blank()
  )

# ---- SAVE PLOT CORRECTLY ----
ggsave(
  filename = paste0(file.rjobs, "/pi_boxplot_3traits_scenario", scenarioIndex, ".png"),
  plot = p,
  width = 10,
  height = 6,
  dpi = 300
)

# ------------------------------------------------------------
# 8) Summary table
# ------------------------------------------------------------
summary_table <- data_long %>%
  group_by(Parameter) %>%
  summarise(
    Mean = mean(Value),
    SD   = sd(Value),
    Median = median(Value),
    Min    = min(Value),
    Max    = max(Value),
    .groups = "drop"
  ) %>%
  left_join(pi_values, by = "Parameter")

write.csv(summary_table,
          file = paste0(file.rjobs, "/summary_table3traits_scenario", scenarioIndex, ".csv"),
          row.names = FALSE)

cat(">>> DONE: Boxplot + Summary table created.\n")
