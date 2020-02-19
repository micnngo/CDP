### CDP: Condensed 20 Newsgroup data 

### Load libraries and functions

source('postDPMManalysis_functions.R')
Rcpp::sourceCpp('reassign_obs_parallel.cpp')


### Load data 
dat <- readRDS('../Data/documents.rds')

# need to remove all empty rows and columns 
which(rowSums(dat) == 0)
which(colSums(dat) == 0)
colnames(dat) <- paste0("D", seq(1:ncol(dat)))
rownames(dat) <- paste0("W", seq(1:nrow(dat)))
dim(dat)

### Run DPMM in Julia 
row_hyper = 1
row_alpha = 10
row_iter = 1000
col_hyper = 1
col_alpha = 100
col_iter = 3000 


# Run on rows (assigns columns)
julia_assign("tdat", t(dat))
julia_assign("row_hyper", row_hyper)
julia_assign("row_alpha", row_alpha)
julia_assign("row_iter", row_iter)
row_dp = julia_eval("DPMMSubClusters.fit(Matrix(tdat), 
                    DPMMSubClusters.multinomial_hyper(fill(row_hyper, size(tdat)[1])), 
                    row_alpha, iters = row_iter, seed = 1234)")

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
                    col_alpha, iters = col_iter, seed = 1234)")

col_LL <- get_LL_and_NClust(col_dp)
col_z <- do.call(rbind, col_dp[[8]])



# Prepare file to save 
fheader = paste0("Results/documents_DPMM_rH_", row_hyper, "_rA_", row_alpha, "_cH_", col_hyper, "_cA_", col_alpha)
h5fname = paste0(fheader, ".h5")
h5createFile(h5fname)
h5createGroup(h5fname, "row")
h5createGroup(h5fname, "col")
h5ls(h5fname)  # check file structure 
h5write(as.matrix(row_LL), h5fname, "row/LL")
h5write(row_z, h5fname, "row/iter_z")
h5write(as.matrix(col_LL), h5fname, "col/LL")
h5write(col_z, h5fname, "col/iter_z")


### Calculate MAP and obtain initial cluster assignments z^r and z^c

# calculate most probable total number of clusters 
top_n = 3
top_r = 2 # choose a number from 1 to top_n
top_c = 1 # choose a number from 1 to top_n
z_r_unsorted = data.frame("row_assignment" = row_z[get_z_idx(row_LL$total_num_clust, top_n = top_n, top = top_r), ])
z_c_unsorted = data.frame("col_assignment" = col_z[get_z_idx(col_LL$total_num_clust, top_n = top_n, top = top_c), ])

rtop = get_top_clusters(row_LL$total_num_clust, top_n = top_n)
ctop = get_top_clusters(col_LL$total_num_clust, top_n = top_n)

# rename z to be in order of ascending labels by convention and then save 
z_r = z_name_convention(z_r_unsorted)
z_c = z_name_convention(z_c_unsorted)


starttime = proc.time()
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

# need to split into correct dimensions
z_r_split <- splitvec(unlist(z_r_updated), rowSums(dat))
z_c_split <- splitvec(unlist(z_c_updated), colSums(dat))

# tabulate updated assignments 
r_table_count = lapply(z_r_split, function(x) tabulate(x, nbins = max(unique(z_r$row_assignment))) )
c_table_count = lapply(z_c_split, function(x) tabulate(x, nbins = max(unique(z_c$col_assignment))) )

# Calculate phi and theta
# assume that alpha = 0 for both phi_r and phi_c
phi_r = as.data.frame(do.call(rbind, lapply(r_table_count, function(x) x/sum(x))))
phi_c = as.data.frame(do.call(rbind, lapply(c_table_count, function(x) x/sum(x))))

z_r_expand = splitz(z_r_updated, dat, 'col')
z_c_expand = splitz(z_c_updated, dat, 'row')

# get frequencies for each pairing and transform into theta 
freq_listoflists = calc_frequency_list(z_r_expand, z_c_expand, dat)
freq_df = get_freq_list(freq_listoflists)
theta_df = get_freq_table(freq_df)
theta_table = theta_into_table(theta_df, "Prob")


endtime = proc.time()
print(endtime - starttime)


### Save intermediary steps 
saveRDS(z_r, paste0(fheader, "_row_z_MAP_", rtop$NumClusters[top_r], ".rds"))
saveRDS(z_c, paste0(fheader, "_col_z_MAP_", ctop$NumClusters[top_c], ".rds"))

saveRDS(z_r_list, paste0(fheader, "_row_z_list_MAP_", rtop$NumClusters[top], ".rds"))
saveRDS(z_c_list, paste0(fheader, "_col_z_list_MAP_", ctop$NumClusters[top], ".rds"))

saveRDS(z_r_updated, paste0(fheader, "_row_z_updated_MAP_", rtop$NumClusters[top], ".rds"))
saveRDS(z_c_updated, paste0(fheader, "_col_z_updated_MAP_", ctop$NumClusters[top], ".rds"))

saveRDS(r_table_count, paste0(fheader, "_row_tabulated_MAP_", rtop$NumClusters[top], ".rds"))
saveRDS(c_table_count, paste0(fheader, "_col_tabulated_MAP_", ctop$NumClusters[top], ".rds"))

saveRDS(phi_r, paste0(fheader, "_row_phi_MAP_", rtop$NumClusters[top], ".rds"))
saveRDS(phi_c, paste0(fheader, "_col_phi_MAP_", ctop$NumClusters[top], ".rds"))

saveRDS(z_r_expand, paste0(fheader, "_row_expanded_MAP_", rtop$NumClusters[top], ".rds"))
saveRDS(z_c_expand, paste0(fheader, "_col_expanded_MAP_", ctop$NumClusters[top], ".rds"))

saveRDS(theta_df, paste0(fheader, "_theta_Kr_", rtop$NumClusters[top_r], "_Kc_", ctop$NumClusters[top_c], ".rds"))


### Visualizations 
pdir = "Plots/documents/"
pheader = paste0("documents_DPMM_rH_", row_hyper, "_rA_", row_alpha, "_cH_", col_hyper, "_cA_", col_alpha, "_")

# plots histogram of total num of clusters and log likelihood 
plot_HGLL(row_LL, pdir, paste0(pheader, "row_"))
plot_HGLL(col_LL, pdir, paste0(pheader, "col_"))

# plots heatmap of theta 
theta_fname = paste0(pdir, pheader, "theta_Kr_", rtop$NumClusters[top_r], "_Kc_", ctop$NumClusters[top_c], "_heatmap.pdf")
plot_theta(theta_df, "Prob", theta_fname)

### Post-Analysis
# top biclusters  
theta_df %>% arrange(desc(Prob))

rgroup = most_probable_z(phi_r)
cgroup = most_probable_z(phi_c)

# example of one bicluster
dat[get_group_idx(rgroup, 1), get_group_idx(cgroup, 2)]
