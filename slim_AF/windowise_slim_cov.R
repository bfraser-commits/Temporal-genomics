library(data.table)
library(dplyr)

# Set constants
min_gen <- 1
mid_gen <- 6
max_gen <- 12
wind_size <- 50000

# Define where your files are
path <- "outputs_vcf_min_outs/"

# Run for 100 iterations
all_winds <- lapply(1:100, function(i) {
  # Read the file
  file_name <- paste0(path, "TY_test_run_iter", i, ".res")
  if (!file.exists(file_name)) return(NULL)  # Skip missing files
  
  test <- data.frame(fread(file_name, fill = TRUE, skip = 15))
  colnames(test) <- c("pos", "freq", "gen")
  
  # Filter generations
  test <- test[test$gen %in% c(min_gen, mid_gen, max_gen),]
  
  # Keep mutations present in all 3 generations
  valid_mut_ids <- names(table(test$pos)[table(test$pos) == 3])
  test <- test[test$pos %in% valid_mut_ids,]
  
  if (nrow(test) == 0) return(NULL)  # Skip empty result
  
  # Reshape to wide format
  test_df <- reshape(test, timevar = "gen", idvar = "pos", direction = "wide")
  
  #remove fixed
  test_df <- test_df[!(test_df$freq.1 == 1 & test_df$freq.12 == 1), ]
  #I_slim <- I_slim[I_slim$p1 >= 0.05 & I_slim$p1 <= 0.95, ]
  
  # Calculate selection coefficients
  test_df$sel_coef_1 <- test_df$freq.6 - test_df$freq.1
  test_df$sel_coef_2 <- test_df$freq.12 - test_df$freq.6
  test_df$sel_coef_total <- test_df$freq.12 - test_df$freq.1
  
  # Define windows
  winds1 <- seq(0, max(test_df$pos), by = wind_size)
  winds2 <- winds1 + wind_size
  
  # Compute summary stats per window
  wind_df <- data.frame(rbindlist(lapply(1:length(winds2), function(y) {
    tmp2 <- test_df[test_df$pos <= winds2[y] & test_df$pos >= winds1[y],]
    
    if (nrow(tmp2) > 0) {
      sel1_tmp <- tmp2$sel_coef_1
      sel2_tmp <- tmp2$sel_coef_2
      n <- length(sel1_tmp)
      x_bar <- mean(sel1_tmp)
      y_bar <- mean(sel2_tmp)
      cov_xy <- sum((sel1_tmp - x_bar) * (sel2_tmp - y_bar)) / (n - 1)
      out <- data.frame(cov = cov_xy, mean_sel1 = x_bar, mean_sel2 = y_bar, site_count = n)
    } else {
      out <- data.frame(cov = NA, mean_sel1 = NA, mean_sel2 = NA, site_count = NA)
    }
    
    out$window <- y
    out$BP1 <- as.integer(winds1[y]) + 1
    out$BP2 <- as.integer(winds2[y])
    out$comp <- 'TA'
    return(out)
  })))
  
  wind_df$iteration <- i  # Add iteration number
  return(wind_df)
})

# Combine all iterations into one data frame
final_df <- bind_rows(all_winds)
write.table(final_df, "TA_vcf_slim_winds_cov.txt", sep = "\t", row.names=F, quote=F)