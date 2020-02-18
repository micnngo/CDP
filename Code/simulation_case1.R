library(dplyr)
library(reticulate)
library(Matrix)
library(ggplot2)
source('simulation_functions.R')
source('fit_CDP.R')
use_python('/path/to/python')
source_python('fit_dt2b.py')
set.seed(123)

# Number of reps of each simulation sestting
n_reps <- 1000

# Simulation parameter values.
N <- 4000
case <- 1
total_cluster_prob <- 0.8


# Data frame to store simulation results.

results <- data.frame(expand.grid(total_cluster_prob = total_cluster_prob, 
                                  case = case, N = N, rep = 1:n_reps),
                      jaccard_plaid = NA, jaccard_bccc = NA, 
                      jaccard_bimax = NA, jaccard_xmotifs = NA,
                      jaccard_spectral = NA, jaccard_dt2b = NA, jaccard_CDP = NA)


start <- proc.time()
for (i in 1:nrow(results)) {
  # Get true bicluster structure.
  df_truth <- get_truth(results$case[i])
  dat <- get_example_data(N = results$N[i], total_cluster_prob = results$total_cluster_prob[i], case_number = results$case[i])
  
  # Compute Jaccard similarity for each biclustering method.
  res_plaid <- tidy_biclust_output(biclust(log(dat + 1), method=BCPlaid(), verbose = FALSE, fit.model = y ~ m + a + b + a:b))
  results$jaccard_plaid[i] <- score(res_plaid, df_truth)
  
  res_bccc <- tidy_biclust_output(biclust(scale(dat), method=BCCC(), number = 5))
  results$jaccard_bccc[i] <- score(res_bccc, df_truth)

  res_bimax <- tidy_biclust_output(biclust(log(dat + 1), method=BCBimax()))
  results$jaccard_bimax[i] <- score(res_bimax, df_truth)
  
  res_xmotifs <- tidy_biclust_output(biclust(dat, method = BCXmotifs()))
  results$jaccard_xmotifs[i] <- score(res_xmotifs, df_truth)
  
  res_spectral <- tidy_biclust_output(biclust(dat, method = BCSpectral()))
  results$jaccard_spectral <- score(res_spectral, df_truth)
  
  res_dt2b <- tidy_dt2b(dat, 5, 5)
  results$jaccard_dt2b[i] <- score(res_dt2b, df_truth)
  
  res_CDP <- fit_CDP(dat, row_alpha = 1, col_alpha = 1,
                           row_hyper = 10, col_hyper = 10, threshold_type = "mean+sd")
  results$jaccard_CDP[i] <- score(df_truth, res_CDP)
}
end <- proc.time()
print(end - start)

# Summarize results
results %>% 
  group_by(case) %>% 
  summarize(mean_plaid = mean(jaccard_plaid), 
            mean_bccc = mean(jaccard_bccc), 
            mean_bimax = mean(jaccard_bimax),
            mean_spectral = mean(jaccard_spectral),
            mean_xmotifs = mean(jaccard_xmotifs),
            mean_dt2b = mean(jaccard_dt2b),
            mean_CDP = mean(jaccard_CDP),
            sd_plaid = sd(jaccard_plaid), 
            sd_bccc = sd(jaccard_bccc), 
            sd_bimax = sd(jaccard_bimax), 
            sd_xmotifs = sd(jaccard_xmotifs),
            sd_spectral = sd(jaccard_spectral),
            sd_dt2b = sd(jaccard_dt2b),
            sd_CDP = sd(jaccard_CDP))
