# Load libraries
lib <- c("data.table","ggplot2","dplyr","parallel","Rfast","cowplot")
lapply(lib, library, character.only = TRUE)

# results directory
sim_res <- "outputs_vcf_min_outs"
fibr_res <- list.files(sim_res)

# Modified selection coefficient calculation for three generations
calc_sel_coef_BF <- function(x){
  # Extract frequencies
  p1 <- x[x$generation == min(x$generation), "freq"]
  pmid <- x[x$generation == sort(unique(x$generation))[2], "freq"]
  pt <- x[x$generation == max(x$generation), "freq"]
  q1 <- 1 - p1
  qmid <- 1 - pmid
  qt <- 1 - pt

  # Time intervals
  gen_order <- sort(unique(x$generation))
  tau_1_mid <- gen_order[2] - gen_order[1]
  tau_mid_t <- gen_order[3] - gen_order[2]

  # Selection coefficients
  s_1_mid <- (2 / tau_1_mid) * log((pmid * q1) / (qmid * p1))
  s_mid_t <- (2 / tau_mid_t) * log((pt * qmid) / (qt * pmid))

  return(data.frame(
    mutation_id = unique(x$mutation_id),
    p1 = p1,
    pmid = pmid,
    pt = pt,
    sel_coef_1_mid = s_1_mid,
    sel_coef_mid_t = s_mid_t
  ))
}

#### Calculate neutral selection coefficients ####
#intro_pops <- c("LL","UL","CA","TY")
intro_pops <- "TY"

if (!file.exists("outputs_test/neutral_sim_13_18.rds")) {
  pop_neutral_SC <- lapply(intro_pops, function(pop) {
    print(paste0("STARTING ", pop))
    pop_runs <- grep(pop, fibr_res, value = TRUE)

    # Set mid generation based on population
    mid_gen <- ifelse(pop %in% c("CA", "TY"), 6, 8)

    run_res <- data.frame(rbindlist(pbmcapply::pbmclapply(pop_runs, function(run) {
      tmp <- data.frame(fread(paste0(sim_res, "/", run), fill = TRUE, skip = 15))
      print(paste0(sim_res, "/", run))
      colnames(tmp) <- c("mutation_id", "freq", "generation")
      tmp$run <- run

      min_gen <- 1
      max_gen <- max(tmp$generation)

      # Subset to the required generations
      tmp <- tmp[tmp$generation %in% c(min_gen, mid_gen, max_gen),]

      # Keep only mutations present in all three generations
      mut_counts <- table(tmp$mutation_id)
      valid_mut_ids <- names(mut_counts[mut_counts == 3])
      tmp <- tmp[tmp$mutation_id %in% valid_mut_ids,]

      # Calculate selection coefficients
      sel_coef_res <- data.frame(rbindlist(lapply(unique(tmp$mutation_id), function(x) {
        calc_sel_coef_BF(tmp[tmp$mutation_id == x, ])
      })))

      return(sel_coef_res)
    }, mc.cores = 4)))

    run_res$pop <- pop
    return(run_res)
  })
  saveRDS(pop_neutral_SC, "outputs_vcf_min/neutral_sim_13_18_TA.rds")
} else {
  pop_neutral_SC <- readRDS("outputs_vcf_min/neutral_sim_13_18_TA.rds")
}

# Combine results
all_pop_sims <- data.frame(rbindlist(pop_neutral_SC))
write.table(all_pop_sims, file = "all_pop_sims_min_outs_13_18_vcf_TA.txt", sep = "\t", row.names = FALSE, quote = FALSE)

