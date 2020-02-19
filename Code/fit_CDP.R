source('postDPMManalysis_functions.R')
Rcpp::sourceCpp('reassign_obs_parallel.cpp')
julia_setup(JULIA_HOME = "/path/to/julia")

fit_CDP <- function(dat, row_alpha = 0.0001, col_alpha = 0.0001, 
                       row_hyper = 0.01, col_hyper = 0.01,
                       threshold_type = "mean") {
  ### Run DPMM in Julia 
  
  julia_library("DPMMSubClusters")
  julia_library("Random")
  julia_library("DataFrames")
  
  totalstarttime = proc.time()
  
#  row_hyper = 0.01
#  row_alpha = 0.0001
  row_iter = 1000
#  col_hyper = 0.01
#  col_alpha = 0.0001
  col_iter = 1000
  
  # Run on rows (assigns columns)
  
  julia_assign("tdat", t(dat))
  julia_assign("row_hyper", row_hyper)
  julia_assign("row_alpha", row_alpha)
  julia_assign("row_iter", row_iter)
  row_dp = julia_eval("DPMMSubClusters.fit(Matrix(tdat), 
                      DPMMSubClusters.multinomial_hyper(fill(row_hyper, size(tdat)[1])), 
                      row_alpha, verbose = false, iters = row_iter, seed = 1234)")
  
  # save log likelihood and total number of clusters 
  row_LL <- get_LL_and_NClust(row_dp)
  
  # save per iteration cluster assignments 
  row_z <- do.call(rbind, row_dp[[8]])
  
  # Run on columns (assigns rows)
  
  julia_assign("dat", dat)
  julia_assign("col_hyper", col_hyper)
  julia_assign("col_alpha", col_alpha)
  julia_assign("col_iter", col_iter)
  col_dp = julia_eval("DPMMSubClusters.fit(Matrix(dat), 
                      DPMMSubClusters.multinomial_hyper(fill(col_hyper, size(dat)[1])), 
                      col_alpha, verbose = false, iters = col_iter, seed = 1234)")
  
  col_LL <- get_LL_and_NClust(col_dp)
  col_z <- do.call(rbind, col_dp[[8]])
  
  ### Calculate MAP and obtain initial cluster assignments z^r and z^c
  
  # calculate most probable total number of clusters 
  top_n = 3
  top = 1 # choose a number from 1 to top_n
  
  z_r_unsorted = data.frame("row_assignment" = row_z[get_z_idx(row_LL$total_num_clust, top_n = top_n, top = top), ])
  z_c_unsorted = data.frame("col_assignment" = col_z[get_z_idx(col_LL$total_num_clust, top_n = top_n, top = top), ])
  get_z_idx(col_LL$total_num_clust, top_n = top_n, top = top)
  rtop = get_top_clusters(row_LL$total_num_clust, top_n = top_n)
  ctop = get_top_clusters(col_LL$total_num_clust, top_n = top_n)
  
  # rename z to be in order of ascending labels by convention and then save 
  z_r = z_name_convention(z_r_unsorted)
  z_c = z_name_convention(z_c_unsorted)
  
  ### Calculate theta, phi_r and phi_c 
  # update assignments using data as weights 
  nCores <- detectCores() - 2
  cl <- makeCluster(nCores)
  clusterExport(cl, varlist=c("z_r", "z_c", "dat"), envir=environment())
  z_r_list = parApply(cl=cl, dat, 2, function(x) rep(z_r$row_assignment, x))
  z_c_list = parApply(cl=cl, dat, 1, function(x) rep(z_c$col_assignment, x))
  stopCluster(cl)
  
  z_r_updated = z_list_reassign(z_r_list, max(unique(z_r$row_assignment)), 100)
  z_c_updated = z_list_reassign(z_c_list, max(unique(z_c$col_assignment)), 100)
  
  # tabulate updated assignments 
  r_table_count = lapply(z_r_updated, function(x) tabulate(x, nbins = max(unique(z_r$row_assignment))) )
  c_table_count = lapply(z_c_updated, function(x) tabulate(x, nbins = max(unique(z_c$col_assignment))) )
  
  # Calculate phi and theta
  phi_r = as.data.frame(do.call(rbind, lapply(r_table_count, function(x) x/sum(x))))
  phi_c = as.data.frame(do.call(rbind, lapply(c_table_count, function(x) x/sum(x))))
  
  z_r_expand = splitz(z_r_updated, dat, 'col')
  z_c_expand = splitz(z_c_updated, dat, 'row')
  
  # get frequencies for each pairing and transform into theta 
  freq_listoflists = calc_frequency_list(z_r_expand, z_c_expand, dat)
  freq_df = get_freq_list(freq_listoflists)
  theta_df = get_freq_table(freq_df)
  theta_table = theta_into_table(theta_df, "Prob")
  
  ### Accuracy
  cluster_list_c <- apply(phi_c, 1, function(x) which(x == max(x)))
  cluster_list_r <- apply(phi_r, 1, function(x) which(x == max(x)))
  
  row_clusters <- data.frame("z_r" = unlist(apply(phi_r, 1, function(x) which(x == max(x)))), 
                             "row" = rep(1:nrow(phi_r), sapply(cluster_list_r, length)))
  col_clusters <- data.frame("z_c" = unlist(apply(phi_c, 1, function(x) which(x == max(x)))),
                             "col" = rep(1:nrow(phi_c), sapply(cluster_list_c, length)))
  
  df_clusters <- data.frame(expand.grid(row = row_clusters[, 2], col = col_clusters[, 2]))
  df_clusters$z_r <- row_clusters[, 1]
  df_clusters$z_c <- rep(col_clusters[, 1], each = nrow(row_clusters))
  df_clusters$bicluster_id <- NA
  
  df_clusters$bicluster_id <- df_clusters$z_c*100 + df_clusters$z_r
  
  for (i in 1:nrow(df_clusters)) {
    df_clusters$count[i] <- dat[df_clusters$row[i], df_clusters$col[i]]
  }
  
  mean_cluster_counts <- df_clusters %>% 
    group_by(bicluster_id) %>% 
    dplyr::summarise(mean_count = mean(count))
  
  overall_mean_count <- mean(mean_cluster_counts$mean_count)
  overall_mean_sd <- sd(mean_cluster_counts$mean_count)
  if (threshold_type == "mean") {
    threshold <- overall_mean_count
  } else if (threshold_type == "mean+sd") {
    threshold <- overall_mean_count + overall_mean_sd  
  } else if (threshold_type == "mean+sd/2") {
    threshold <- overall_mean_count + overall_mean_sd/2
  }
  
  for (i in 1:nrow(df_clusters)) {
    mean_cluster_count <- mean_cluster_counts$mean_count[mean_cluster_counts$bicluster_id == df_clusters$bicluster_id[i]]
    if (mean_cluster_count < threshold) {
      df_clusters$bicluster_id[i] <- 0
    }
  }
  df_clusters <- df_clusters %>% filter(bicluster_id != 0)
  
  df_results <- tibble("gene" = paste0("R", as.character(df_clusters$row)), "cell" = paste0("C", df_clusters$col),
                       "bicluster_id" = df_clusters$bicluster_id)
  return(df_results)
}
