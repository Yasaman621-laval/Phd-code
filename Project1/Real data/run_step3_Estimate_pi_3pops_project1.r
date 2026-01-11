rm(list=ls())
library(SummaryLasso)


# Constants and Directory Setup
testset <- 2
savename <- "new_densityU_3pops_project1"
penalty <- "RealmixLOG"

  dirOutput = "/home/yatah3/projects/def-thchlava/yatah3/real-data/project1/"

# Population and Trait Settings
popvec <- c("AFR", "EAS", "EUR")

rho_vec <- c(seq(0, 0.9, 0.1), 0.95)

# Load tau values
  Taufile =paste0(dirOutput,popvec[1],"_",popvec[2], "-", popvec[3],"Tau_info_v2.RData")
  load(file=Taufile)

K <- 3
epsilon <- 1e-8
  
  output_sub_folder_uDensity <- paste0(dirOutput, savename, "/")
  
  output_sub_folder_Beta <- paste0(dirOutput, "/Allchr_beta/")
  dir.create(output_sub_folder_Beta, recursive = TRUE, showWarnings = FALSE)
  
  largest_values <- numeric()
  largest_indices <- numeric()
  largest_rhos <- numeric()
  largest_tauuse <- numeric()
  
  for (tt in 1:length(AbsTauvec)) {
    tauuse <- AbsTauvec[tt]
    print(paste0("Processing tauuse=", tauuse))
    Allrho_ZMatrix_list <- list()
    
    for (rr in 1:length(rho_vec)) {
      rho <- rho_vec[rr]
      print(paste0("rho=", rho))
      
      AllBeta <- NULL
      
      for (chr in 1:22) {
        print(paste0("chr=", chr))
        savefile1 <- paste0(output_sub_folder_uDensity, penalty, "chr", chr, 
                             "usedtrait_1,2,3", "warmStart1", 
                             "Zscale1", "singleStart1", 
                             "tauuse", tauuse, "_DensityU.RData")
        
        if (file.exists(savefile1) && file.info(savefile1)$size > 0) {
          load(savefile1)
        } else {
          print(paste0("Missing file for chr=", chr, ", tauuse=", tauuse))
          next
        }
        
if (length(Allrho_WBMatrix_list) < rr) {print(paste0("Allrho_WBMatrix_list does not contain enough elements, rr=", rr)); next}

        if (chr == 1) {
          AllBeta <- t(Allrho_WBMatrix_list[[rr]])
        } else {
          AllBeta <- cbind(AllBeta, t(Allrho_WBMatrix_list[[rr]]))
        }
        
        rm(Allrho_WBMatrix_list)
        gc()
      }
      
      savefile <- paste0(output_sub_folder_Beta, "usedtrait_1,2,3","allchrs", ".RData")
      Allrho_ZMatrix_list[[rr]] <- AllBeta
      rm(AllBeta)
      gc()
    }
    
    save(Allrho_ZMatrix_list, file = savefile)
    rm(Allrho_ZMatrix_list)
    gc()
    
    results_list <- list()
    
    for (rr in 1:length(rho_vec)) {
      rho <- rho_vec[rr]
      
      if (!file.exists(savefile)) {
        print(paste0("Missing tauuse=", tauuse))
        next
      }
      
      load(savefile)
      combined_matrix <- Allrho_ZMatrix_list[[rr]]
      adjusted_matrix <- ifelse(combined_matrix == 0, combined_matrix + epsilon, combined_matrix)
      sum_of_logs <- rowSums(log(adjusted_matrix))
      
      max_index <- which.max(sum_of_logs)
      max_sum_of_logs <- sum_of_logs[max_index]
      
      results_list[[rr]] <- list(
        Row_index = max_index,
        Largest_sum_of_logs = max_sum_of_logs,
        Corresponding_rho = rho,
        Corresponding_tauuse = tauuse
      )
    }
    
    largest_values_tt <- sapply(results_list, function(result) result$Largest_sum_of_logs)
    largest_indices_tt <- sapply(results_list, function(result) result$Row_index)
    largest_rhos_tt <- sapply(results_list, function(result) result$Corresponding_rho)
    largest_tauuse_tt <- sapply(results_list, function(result) result$Corresponding_tauuse)
    
    max_largest_value_tt <- max(largest_values_tt)
    max_index_tt <- which(largest_values_tt == max_largest_value_tt)
    corresponding_rho_tt <- largest_rhos_tt[max_index_tt]
    corresponding_tauuse <- largest_tauuse_tt[max_index_tt]
    
    largest_values <- c(largest_values, max_largest_value_tt)
    largest_indices <- c(largest_indices, largest_indices_tt[max_index_tt])
    largest_rhos <- c(largest_rhos, corresponding_rho_tt)
    largest_tauuse <- c(largest_tauuse, corresponding_tauuse)
  }
  
  overall_max_index <- which.max(largest_values)
  overall_largest_value <- largest_values[overall_max_index]
  overall_largest_index <- largest_indices[overall_max_index]
  overall_corresponding_rho <- largest_rhos[overall_max_index]
  overall_tauuse <- largest_tauuse[overall_max_index]
  
 data.frame(
    Largest_Sum_Of_Logs = overall_largest_value,
    Largest_Index = overall_largest_index,
    Corresponding_Rho = overall_corresponding_rho,
    Corresponding_Tauuse = overall_tauuse
  )

