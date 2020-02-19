##########
### Post-DPMM Analysis Functions
##########

### Libraries 

package_list <- c("ggplot2", "gridExtra", "JuliaCall", "parallel", 
                  "plyr", "dplyr", "reshape2", "Rcpp", "RcppParallel")
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bio_package_list <- c("rhdf5")
new_bio_packages <- bio_package_list[!(bio_package_list %in% installed.packages()[,"Package"])]
if(length(new_bio_packages)) BiocManager::install(new_bio_packages)

library(ggplot2)
library(gridExtra)
library(JuliaCall)
library(rhdf5)    # for saving Julia output
library(parallel)
library(plyr)
library(dplyr)
library(reshape2) # for acast
library(Rcpp)
library(RcppParallel)


julia_library("DPMMSubClusters")
julia_library("Random")
julia_library("DataFrames")
julia_library("HDF5")


### Gets LogLikelihood and Total Num of Clusters per Iteration 
get_LL_and_NClust <- function(DPMMfitresult){
  # input: JuliaTuple from DPMMSubClusters.fit function (NOT dp_parallel function)
  # gets the Log Likelihood and Total Number of Clusters per iteration
  # returns data frame 
  LL <- unlist(DPMMfitresult[[6]])
  NC <- unlist(DPMMfitresult[[7]])
  df <- data.frame("loglikelihood" = LL, "total_num_clust" = NC)
  return(df)
}

### Plots Histogram of Total Num of Clusters and LogLikelihood per Iteration 
plot_HGLL <- function(df, fdir, fheader){
  # "Plots/DPDT2B/Marques_row_Post_NumClust.pdf"
  # fdir: "Plots/DPDT2B/"
  # fheader: "Marques_row_"
  fontsize = 25
  HG = ggplot(df, aes(x=total_num_clust)) + 
    geom_histogram(color="black", fill="blue", alpha=0.7, binwidth=1) + 
    #scale_x_continuous(name = "Number of Clusters", limits=c(0, 2.1)) +
    labs(y = "Count", x = "Number of Clusters") + 
    theme(axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize), text = element_text(size=fontsize), axis.title = element_text(size=fontsize))
  
  LL = ggplot(data=df, aes(x=as.numeric(rownames(df)), y=loglikelihood, group=1)) +
    geom_path(color="black")  + 
    labs(y = "Log-likelihood", x = "Iteration") + 
    theme(axis.text.x = element_text(size=fontsize), axis.text.y = element_text(size=fontsize), text = element_text(size=fontsize), axis.title = element_text(size=fontsize))
  
  plts = grid.arrange(HG, LL, ncol=1, nrow=2)
  ggsave(plts, filename = paste0(fdir, fheader, "Post_LL_NumClust.pdf"), scale=1.5, dpi = 'retina', width = 11.5, height = 8, units="in", device="pdf")  
  ggsave(HG, filename = paste0(fdir, fheader, "Post_NumClust.pdf"), scale=1.5, dpi = 'retina', width = 11.5, height = 8, units="in", device="pdf")  
  ggsave(LL, filename = paste0(fdir, fheader, "Post_LL.pdf"), scale=1.5, dpi = 'retina', width = 11.5, height = 8, units="in", device="pdf")  
  return(plts)
}



### Gets the top_n most probable total number of clusters
get_top_clusters <- function(NC, top_n=3){
  # input
  #   NC : vector of total number of clusters at each iteration
  #   top_n : returns the n most probable number of clusters 
  # returns (top_n x 2) dataframe showing the max number of clusters (K) and its probability 
  tbl = table(NC)
  tbl_ordered = tbl[order(-tbl)]/sum(tbl)
  df = data.frame(tbl_ordered[1:top_n])
  
  if(ncol(df) == 1){
    colnames(df) = "NumClusters"
    df$Prob = 1 
  } else{
    colnames(df) = c("NumClusters", "Prob")
  }
  return(df)
}

### Gets the last index of the DPMM iteration where Total Num of Clusters = top 
get_z_idx <- function(NC, top_n, top){
  # input
  #   NC: total number of clusters per iteration
  #   top: int from 1 to top_n; rtop[top] chooses K for rows 
  most_probable_MAP = get_top_clusters(NC, top_n)
  print(most_probable_MAP)
  indexes = which(NC == most_probable_MAP$NumClusters[top])
  lastidx = indexes[length(indexes)]
  return(lastidx)
}


### Renames z_r/z_c labels to be in order (by convention z = 1, ..., K)
z_name_convention <- function(z_updated){
  z_updated = as.matrix(z_updated, nrow=dim(z_updated)[1], ncol=dim(z_updated)[2])
  
  # set all cluster labels to one if all obs in same cluster 
  if(length(unique(z_updated)) == 1){
    z_new <- matrix(rep(1, length(z_updated)), nrow=dim(z_updated)[1], ncol=dim(z_updated)[2])
    colnames(z_new) <- colnames(z_updated)
  } else{
    # otherwise reassign cluster labels so they go from 1, 2, ... K
    K_max = max(unique(z_updated))
    original_order = unique(z_updated)
    z_new <- plyr::mapvalues(z_updated, original_order, seq(1, K_max))
  }
  
  return(as.data.frame(z_new))
}


##### 
### Splits lists to calculate phi_r, phi_c and theta 
#####

# splits each expanded list of assignments into sublists of length according to data matrix 
splitvec <- function(vec, splitlen){
  sv = split(vec, rep(seq_along(splitlen), splitlen))
  return(sv)
}

# apply the split command to either the row or col of list of assignments 
splitz <- function(zlist, df, zdir='row'){
  dir = tolower(zdir)
  if(dir == "row"){
    sz = lapply(1:nrow(df), function(x) splitvec(zlist[[x]], df[x,]))
    names(sz) = rownames(df)
  } else if (dir == "col"){
    sz = lapply(1:ncol(df), function(x) splitvec(zlist[[x]], df[,x]))
    names(sz) = colnames(df)
  } else
    print('not valid direction')
  return(sz)
}


### Calculates the frequency of each row and col cluster assignment 
calc_frequency_list <- function(z_r_exp, z_c_exp, df){
  frequency_list <- list()
  for(i in 1:length(z_c_exp)){
    # apply tally function
    #print(i)
    frequency_list[[i]] <- applytally2(z_c_exp[[i]], z_r_exp, df, i)
  }
  return(frequency_list)
}

# creates a list of tallied counts (not condensed)
applytally2 <- function(sublist, other_list, df, current_idx){
  # the other list MUST BE z_r_expand 
  sublist_length = length(sublist)
  other_list_names = colnames(df)[as.integer(names(sublist))]
  
  tally_list <- list()
  for(i in 1:sublist_length){
    #print(length(sublist[[i]]))
    #print(length(other_list[[other_list_names[i]]][[current_idx]]))
    #print(i)
    #print(other_list_names[i])
    tally_list[[i]] <- match_columns(sublist[[i]], other_list[[other_list_names[i]]], current_idx)
    #print(tally_list)
  }
  return(tally_list)
}

# tabulates the row and col vectors 
match_columns <- function(ref_vector, other_sublist, current_idx){
  # tables pairings between two vectors
  other_vector = other_sublist[[as.character(current_idx)]]
  #print(length(other_vector))
  #print(head(data.frame(rvec=other_vector, cvec=ref_vector), 10))
  #print(table(data.frame(rvec=other_vector, cvec=ref_vector)))
  tally = data.frame(table(data.frame(rvec=other_vector, cvec=ref_vector)))
  return(tally)
}



### converts nested list of frequencies into a data frame with duplicated entries 
get_freq_list <- function(freq_nestedlists){
  # changes output of lapply where function is applytally to long data frame of freq
  freq_list = lapply(freq_nestedlists, function(x) reshape2::melt(x, id=c("rvec", "cvec")))
  freq_longdf = reshape2::melt(freq_list, id=c("rvec", "cvec", "variable"))
  freq_longdf = freq_longdf[,-5] # drop the L1 column 
  freq_longdf = freq_longdf %>% arrange(rvec, cvec)
  return(freq_longdf)
}

### transforms data frame from get_freq_list into data frame version of theta 
get_freq_table <- function(freq_df){
  # data should be in format: rvec, cvec, variable (Freq), value
  # gives theta in data.frame format 
  freq_table = xtabs(value~. , data = freq_df)
  freq_df = data.frame(freq_table)
  freq_df$Prob = freq_df$Freq/sum(freq_df$Freq)
  freq_df = freq_df[, -3] # drop the variable column  
  return(freq_df)
}

### transforms theta (data frame) into contingency table version 
theta_into_table <- function(theta, valtype = "Freq"){
  # value type: "Freq" or "Prob" 
  # theta format: rvec, cvec, Freq, Prob
  c_table = acast(theta, rvec~cvec, value.var = valtype)
  roworder = as.character(sort(as.integer(rownames(c_table)), decreasing=FALSE))
  colorder = as.character(sort(as.integer(colnames(c_table)), decreasing=FALSE))
  contin_table = c_table[roworder, colorder]
  return(contin_table)
}

### plot heatmap of theta 
plot_theta <- function(theta_df, filltype = "Prob", fname){
  thetafig <- ggplot(theta_df, aes(y=rvec, x=cvec )) + 
    geom_tile(aes_string(fill = filltype), color="grey") +
    scale_x_discrete(limits = as.character(sort(as.integer(levels(theta_df$cvec))))) + 
    scale_y_discrete(limits = as.character(sort(as.integer(levels(theta_df$rvec))))) + 
    labs(y = "Latent Row Var", x = "Latent Col Var", fill = filltype) + 
    theme(axis.text.x = element_text(size=40), axis.text.y = element_text(size=40), aspect.ratio = 1, text = element_text(size=40), axis.title = element_text(size=40), legend.text=element_text(size=40),
          legend.title=element_text(size=40), legend.key.size = unit(3,"line"))
  ggsave(thetafig, filename = fname, scale=2, dpi = 'retina', width = 11, height = 8, units="in", device="pdf")  
  return(thetafig)
}


### determines the most probable cluster assignments
most_probable_z <- function(phi){
  # input: phi_r or phi_c data frame 
  probable_z = unlist(apply(phi, 1, function(x) which(x == max(x)) ))

  return(probable_z)
}


### gets indicies of observations in cluster K = group_number
get_group_idx <- function(unlisted_group, group_number){
  # inputs:
  #   unlisted_group: unlisted array of most probable cluster assignments for each observation (need phi_r/phi_c)
  #   group_number: which cluster is of interest 
  grp = names(which(unlisted_group == group_number))
  idx =  as.numeric( unique(unlist(lapply(sapply(grp, function(x) strsplit(x, ".V")), "[[", 1))) )
  
  return(idx)
}
