library(biclust)

# Tidy the output from biclust(). 
# Returns a df with columns bicluster_id, gene (row), cell (col)
# corresponding to the biclustering given in res.
tidy_biclust_output <- function(res) {
  n_biclusters <- res@Number
  row_counts <- colSums(res@RowxNumber)
  if (n_biclusters == 0) {
    return(NULL)
  } else if (n_biclusters == 1) {
    col_counts <- sum(res@NumberxCol)
  } else {
    col_counts <- rowSums(res@NumberxCol)
  }
  n_bicluster_elements <- row_counts * col_counts
  df_biclusters <- data.frame(bicluster_id = rep(1:n_biclusters, n_bicluster_elements),
                              gene = NA, cell = NA)
  if (n_biclusters > 1) {
    for (i in 1:n_biclusters) {
      df_biclusters[df_biclusters$bicluster_id == i, c("gene", "cell")] <- 
        expand.grid(which(res@RowxNumber[, i]), which(res@NumberxCol[i, ]))
    }
  } else {
    df_biclusters[, c("gene", "cell")] <- expand.grid(which(res@RowxNumber), which(res@NumberxCol))
  }
  df_biclusters$gene <- paste0("R", df_biclusters$gene)
  df_biclusters$cell <- paste0("C", df_biclusters$cell)
  return(df_biclusters)
}

# Return DT2B bicluster results in tidy format.
tidy_dt2b <- function(dat, Kr, Kc) {
  dimnames(dat) <- list(paste0("R", 1:nrow(dat)), paste0("C", 1:ncol(dat)))
  dat_table <- matrix_to_table(dat)
  readr::write_csv('temp.csv', x = dat_table, col_names = FALSE)
  res_dt2b <- fit_dt2b('temp.csv', Kr, Kc)
  res_dt2b <- threshold_clusters(res_dt2b, 
                                gamma_joint = mean(res_dt2b$joint) + sd(res_dt2b$joint),
                                gamma_row = mean(res_dt2b$phi_r) + sd(res_dt2b$phi_r),
                                gamma_col = mean(res_dt2b$phi_c) + sd(res_dt2b$phi_c))
  res_dt2b <- add_labels(res_dt2b)
  df_biclusters <- tidy_biclusters(res_dt2b)
  df_biclusters <- add_labels_tidybiclusters(res_dt2b, df_biclusters) %>% 
    dplyr::select(c("bicluster_id", "cell", "gene"))
  return(df_biclusters)
}

# Generate count data matrix for N total records, R rows, C columns
# and R x C probability matrix prob (with sum(prob) = 1)
gen_data <- function(N, R, C, prob) {
  dat_vec <- rmultinom(1, N, prob = c(prob))
  dat <- matrix(dat_vec, nrow = R, ncol = C)
  return(dat)
}

# Convert count data matrix to longform table
matrix_to_table <- function(dat) {
  d1 <- expand.grid(dimnames(dat))
  d2 <- d1[rep(1:nrow(d1), c(dat)), ]
  dat_table <- d2[order(d2$Var1), ]
  row.names(dat_table) <- NULL
  return(dat_table)
}

# Select a bicluster with given id from a tidy biclustering.
select_bicluster <- function(id, df) {
  df %>% dplyr::filter(bicluster_id == id) %>% dplyr::select(cell, gene)
}

# Compute Jaccard score between two biclusters.
jaccard <- function(bicluster1, bicluster2) {
  suppressMessages(nrow(inner_join(bicluster1, bicluster2)) / nrow(full_join(bicluster1, bicluster2)))
}

# Compute max Jaccard score between a given bicluster and a biclustering.
max_jaccard <- function(bicluster, df) {
  current_max <- 0
  unique_ids <- unique(df$bicluster_id)
  for (i in 1:length(unique_ids)) {
    bicluster2 <- select_bicluster(unique_ids[i], df)
    val <- jaccard(bicluster, bicluster2)
    current_max <- ifelse(val < current_max, current_max, val)
  }
  return(current_max)
}

# Directed Jaccard score.
directed_score <- function(df1, df2) {
  result <- 0
  unique_ids <- unique(df1$bicluster_id)
  for (i in 1:length(unique_ids)) {
    b1 <- select_bicluster(unique_ids[i], df1)
    result <- result + max_jaccard(b1, df2)
  }
  result <- result / length(unique(df1$bicluster_id))
  return(result)
}

# Compute Jaccard similarity score between two tidy biclusterings.
score <- function(df1, df2) {
  if (length(df1) == 0 | length(df2) == 0) {
    return(0)
  } else {
    return(min(directed_score(df1, df2), directed_score(df2, df1)))
  }
}

# Return bicluster results in tidy format.
tidy_biclusters <- function(results) {
  results_summary <- as.data.frame(which(results$joint > 0, arr.ind = TRUE))
  topics <- list()
  for (i in 1:nrow(results_summary)) {
    topics[[i]] <- list(row = which(results$phi_r[results_summary$row[i], ] > 0), col = which(results$phi_c[results_summary$col[i], ] > 0))
  }
  
  for (i in 1:nrow(results_summary)) {
    results_summary$n_elements[i] <- length(topics[[i]]$row) * length(topics[[i]]$col)
  }
  n_elements <- sum(results_summary$n_elements)
  n_biclusters <- nrow(results_summary)
  
  df_biclusters <- data.frame(bicluster_id = rep(NA, n_elements), row_topic = NA, col_topic = NA, row = NA, col = NA)
  df_biclusters$bicluster_id <- rep(1:n_biclusters, results_summary$n_elements)
  df_biclusters$row_topic <- results_summary$row[df_biclusters$bicluster_id]
  df_biclusters$col_topic <- results_summary$col[df_biclusters$bicluster_id]
  rowcol_ids <- matrix(nrow = 0, ncol = 2)
  for (i in 1:n_biclusters) {
    rowcol_ids <- rbind(rowcol_ids, expand.grid(topics[[i]]$row, topics[[i]]$col))
  }
  df_biclusters$row <- rowcol_ids[, 1]
  df_biclusters$col <- rowcol_ids[, 2]
  return(df_biclusters)
}

# Threshold results for DT2B.
threshold_clusters <- function(dt2b_results, gamma_joint, gamma_row, gamma_col) {
  dt2b_results$joint[dt2b_results$joint < gamma_joint] <- 0
  dt2b_results$phi_r[dt2b_results$phi_r < gamma_row] <- 0
  dt2b_results$phi_c[dt2b_results$phi_c < gamma_col] <- 0
  return(dt2b_results)
}

# Add row/col labels for raw DT2B output.
add_labels <- function(results){
  nGenes = ncol(results$phi_r)
  nCells = ncol(results$phi_c)

  genes_idx2val = gsub('[\"]',"", unlist(results$idx2vals)[1:nGenes])
  cells_idx2val = gsub('[\"]',"", unlist(results$idx2vals)[nGenes+1:nCells])
  
  colnames(results$phi_r) <- genes_idx2val
  colnames(results$phi_c) <- cells_idx2val
  
  return(results)
}

# Add row/col labels for tidied DT2B output.
add_labels_tidybiclusters <- function(results, tidy_biclust_df){
  nGenes = ncol(results$phi_r)
  nCells = ncol(results$phi_c)
  
  genes_idx2val = gsub('[\"]',"", unlist(results$idx2vals)[1:nGenes])
  cells_idx2val = gsub('[\"]',"", unlist(results$idx2vals)[nGenes+1:nCells])
  
  tidy_biclust_df$gene = genes_idx2val[tidy_biclust_df$row]
  #tidy_biclust_df$gene = sapply(tidy_biclust_df$gene, toupper) 
  tidy_biclust_df$cell = cells_idx2val[tidy_biclust_df$col]
  
  return(tidy_biclust_df)
}

# Return true bicluster structure in tidy format for specified case.
get_truth <- function(case_number) {
  if (case_number == 1) {
    df_truth <- data.frame(bicluster_id = rep(1, 400), gene = NA, cell = NA)
    df_truth[, c("gene", "cell")] <- rbind(expand.grid(1:20, 1:20))
  } else if (case_number == 2) {
    df_truth <- data.frame(bicluster_id = rep(1:2, c(16, 128), gene = NA, cell = NA))
    df_truth[, c("gene", "cell")] <- rbind(expand.grid(1:4, 1:4), 
                                           expand.grid(5:12, 5:20))
  } else if (case_number == 3) {
    df_truth <- data.frame(bicluster_id = rep(1:4, c(196, 36, 80, 36), gene = NA, cell = NA))
    df_truth[, c("gene", "cell")] <- rbind(expand.grid(1:14, 1:14), 
                                           expand.grid(15:20, 15:20),
                                           expand.grid(21:36, 21:25),
                                           expand.grid(40:50, 32:46))
  } else if (case_number == 4) {
    df_truth <- data.frame(bicluster_id = rep(1:5, c(169, 49, 36, 25, 64), gene = NA, cell = NA))
    df_truth[, c("gene", "cell")] <- rbind(expand.grid(1:13, 1:13), 
                                           expand.grid(14:20, 14:20),
                                           expand.grid(30:35, 30:35),
                                           expand.grid(36:40, 36:40),
                                           expand.grid(43:50, 43:50))
  }
  df_truth$gene <- paste0("R", df_truth$gene)
  df_truth$cell <- paste0("C", df_truth$cell)
  return(df_truth)
}

# Generates count data matrix for each exmample case
get_example_data <- function(N, total_cluster_prob, case_number) {
  if (case_number == 1) {
    R <- 50
    C <- 50
    cluster_indicator <- matrix(0, nrow = R, ncol = C)
    cluster_indicator[1:20, 1:20] <- 1
    cluster_prob <- total_cluster_prob / sum(cluster_indicator)
    non_cluster_prob <- (1 - total_cluster_prob) / sum(cluster_indicator == 0)
    prob <- cluster_indicator * cluster_prob
    prob[prob == 0] <- non_cluster_prob
    dat <- gen_data(N, R, C, prob)
  } else if (case_number == 2) {
    R <- 20
    C <- 20
    cluster_indicator <- matrix(0, nrow = R, ncol = C)
    cluster_indicator[1:4, 1:4] <- 1
    cluster_indicator[5:12, 5:20] <- 1
    cluster_prob <- total_cluster_prob / sum(cluster_indicator)
    non_cluster_prob <- (1 - total_cluster_prob) / sum(cluster_indicator == 0)
    prob <- cluster_indicator * cluster_prob
    prob[prob == 0] <- non_cluster_prob
    dat <- gen_data(N, R, C, prob)
  } else if (case_number == 3) {
    R <- 50
    C <- 50
    cluster_indicator <- matrix(0, nrow = R, ncol = C)
    cluster_indicator[1:14, 1:14] <- 1
    cluster_indicator[15:20, 15:20] <- 1
    cluster_indicator[21:36, 21:25] <- 1
    cluster_indicator[40:50, 32:46] <- 1

    cluster_prob <- total_cluster_prob / sum(cluster_indicator)
    non_cluster_prob <- (1 - total_cluster_prob) / sum(cluster_indicator == 0)
    prob <- cluster_indicator * cluster_prob
    prob[cluster_indicator == 0] <- non_cluster_prob
    dat <- gen_data(N, R, C, prob)
  } else if (case_number == 4) {
    R <- 100
    C <- 100
    cluster_indicator <- matrix(0, nrow = R, ncol = C)
    cluster_indicator[1:13, 1:13] <- 1
    cluster_indicator[14:20, 14:20] <- 1
    cluster_indicator[30:35, 30:35] <- 1
    cluster_indicator[36:40, 36:40] <- 1
    cluster_indicator[43:50, 43:50] <- 1
    cluster_prob <- total_cluster_prob / sum(cluster_indicator)
    non_cluster_prob <- (1 - total_cluster_prob) / sum(cluster_indicator == 0)
    prob <- cluster_indicator * cluster_prob
    prob[cluster_indicator == 0] <- non_cluster_prob
    dat <- gen_data(N, R, C, prob)
  }
  rownames(dat) <- paste0("R", 1:nrow(dat))
  colnames(dat) <- paste0("C", 1:ncol(dat))
  return(dat)
}
