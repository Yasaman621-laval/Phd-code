rm(list=ls())
library(stringr)


#__________________________________________________________________________-
#WB1
#_______________________________________
# ----------------------------
# 1. List all .RData files
# ----------------------------
r2_dir <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/step3_AUC_computation/R2output"

# Get only WB2 files for usedtrait2,3
files <- list.files(r2_dir, pattern = "usedtrait1,3.*_wB1AUC\\.RData$", full.names = TRUE)

# ----------------------------
# 2. Extract tau and rho from filenames
# ----------------------------

params <- data.frame(
  file = files,
  tau  = as.numeric(str_extract(files, "(?<=tauuse)\\d+\\.?\\d*(e-?\\d+)?")),
  rho  = as.numeric(str_extract(files, "(?<=rho)\\d*\\.?\\d*")),
  stringsAsFactors = FALSE
)


results <- data.frame()

for (i in seq_len(nrow(params))) {
  f <- params$file[i]
  load(f)
  
  # Skip if missing or invalid
  if (!exists("AUCsummary")) next
  auc_vals <- as.numeric(AUCsummary["AUC", ])
  if (all(is.na(auc_vals))) next
  
  # Find best model within this file
k_best <- which.max(auc_vals)
  auc_best <- max(auc_vals, na.rm = TRUE)
  var_best <- AUCsummary["Var", k_best]
  
  results <- rbind(
    results,
    data.frame(
      tau = params$tau[i],
      rho = params$rho[i],
      best_model = k_best,
      AUC = auc_best,
      Var = var_best,
      file = basename(f)
    )
  )
}

# --- Find overall best AUC ---
best_row <- results[which.max(results$AUC), ]

cat("í ¼í¿† Max AUC across all files:\n")
print(best_row)




#__________________________________________________________________________-
#WB2
#_______________________________________
# ----------------------------
# 1. List all .RData files
# ----------------------------
r2_dir <- "/home/yatah3/projects/def-thchlava/yatah3/real-data/step3_AUC_computation/R2output"

# Get only WB2 files for usedtrait2,3
files <- list.files(r2_dir, pattern = "usedtrait1,3.*_wB2AUC\\.RData$", full.names = TRUE)

# ----------------------------
# 2. Extract tau and rho from filenames
# ----------------------------

params <- data.frame(
  file = files,
  tau  = as.numeric(str_extract(files, "(?<=tauuse)\\d+\\.?\\d*(e-?\\d+)?")),
  rho  = as.numeric(str_extract(files, "(?<=rho)\\d*\\.?\\d*")),
  stringsAsFactors = FALSE
)


results <- data.frame()

for (i in seq_len(nrow(params))) {
  f <- params$file[i]
  load(f)
  
  # Skip if missing or invalid
  if (!exists("AUCsummary")) next
  auc_vals <- as.numeric(AUCsummary["AUC", ])
  if (all(is.na(auc_vals))) next
  
  # Find best model within this file
k_best <- which.max(auc_vals)
  auc_best <- max(auc_vals, na.rm = TRUE)
  var_best <- AUCsummary["Var", k_best]
  
  results <- rbind(
    results,
    data.frame(
      tau = params$tau[i],
      rho = params$rho[i],
      best_model = k_best,
      AUC = auc_best,
      Var = var_best,
      file = basename(f)
    )
  )
}

# --- Find overall best AUC ---
best_row <- results[which.max(results$AUC), ]

cat("í ¼í¿† Max AUC across all files:\n")
print(best_row)



