#include "reassign_obs_parallel.h"


// UPDATES ASSIGNMENTS USING DATA AS WEIGHTS 

// Tabulate function, tabulates assignment vector and returns probability of assignment for obs i 

// [[Rcpp::export]]
std::vector<double> tabulate2(const std::vector<double>& x, const unsigned max) {
  std::vector<double> counts(max, 0);
  std::size_t n = x.size();
  for (std::size_t i = 0; i < n; i++) {
    if (x[i] > 0 && x[i] <= max) {
      counts[x[i] - 1]++;
    }
  }
  // bind2nd takes parameter provided and forces it to be the second parameter in the function
  // transform loops through given array and takes value at index/max-1
  std::transform(counts.begin(), counts.end(), counts.begin(), std::bind2nd(std::divides<double>(), max-1));
  return counts;
}

// [[Rcpp::export]]
double c_sample(std::vector<double> p) {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::discrete_distribution<> d(std::begin(p), std::end(p));
  
  double return_val = d(gen);
  
  return return_val + 1;
}

//To call it
// [[Rcpp::export]]
Rcpp::NumericVector testapply(Rcpp::NumericVector x, int numbins) {
  
  std::vector<double> x_new = Rcpp::as<std::vector<double>>(x);
  
  // Make another copy as the representation of z_new
  Rcpp::NumericVector z_new = Rcpp::clone(x);
  
  std::vector<double> prob = tabulate2(x_new, numbins);
  // create the worker
  test_apply3 test_apply3(x_new, numbins, z_new, prob);
  
  // Call it
  parallelFor(0, x_new.size(), test_apply3);
  
  return z_new;
  
}



// Iteratively reassigns each observation using probabilities from calc_prob2

// [[Rcpp::export]]
Rcpp::NumericVector reassign_obs(Rcpp::NumericVector x, int numbins, int niter) {
  std::vector<double> x_new = Rcpp::as<std::vector<double>>(x);
  Rcpp::NumericVector z_new; 
  for (int i=1; i <= niter; i++){
    int N = x.length();
    
    if(N==0){ // if gene or cell has no expression 
      z_new = 0;
    } 
    else if(N==1){
      Rcpp::IntegerVector idx = Rcpp::seq_len(numbins);
      std::vector<double> p = tabulate2(x_new, numbins);
      z_new = c_sample(p);
    }
    else{
      z_new = testapply(x, numbins);
    }
  }
  return z_new;
}

// Loops over list of row/col assignments (z_r_list or z_c_list) and reassigns obs

// [[Rcpp::export]]
Rcpp::List z_list_reassign(Rcpp::List zlist, int numbins, int niter){
  int n = zlist.length();
  Rcpp::List z_updated(n);
  
  for(int i=0; i < n; i++){
    Rcpp::NumericVector v = zlist[i];
    z_updated[i] = reassign_obs(v, numbins, niter);
  }
  return z_updated;
}