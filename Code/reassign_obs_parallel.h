#ifndef REASSIGN_OBS_PARALLEL_H
#define REASSIGN_OBS_PARALLEL_H

#include <Rcpp.h>
#include <algorithm>
#include <RcppParallel.h>
#include <random>
#include <iterator>



std::vector<double> tabulate2(const std::vector<double>& x, const unsigned max);
double c_sample(std::vector<double> p);
Rcpp::NumericVector testapply(Rcpp::NumericVector x, int numbins);

// [[Rcpp::depends(RcppParallel)]]
struct test_apply3 : public RcppParallel::Worker {
  
  // Input NumericalVector
  const std::vector<double> input;
  
  // tabulated stuff
  std::vector<double> tab;
  
  // length of x vector
  int inputlength;
  
  // Output
  RcppParallel::RVector<double> output;
  
  // Constructor
  test_apply3(std::vector<double> input, int inputlength, Rcpp::NumericVector output, std::vector<double> tab)
  : input(input), inputlength(inputlength), output(output), tab(tab) {}
  
  // do the stuff
  void operator() (std::size_t begin, std::size_t end) {
    
    for (std::size_t i = begin; i < end; i++) {
      std::vector<double> p_copy = tab;
      p_copy[input[i]-1] -= 1.0/(inputlength-1);
      output[i] = c_sample(p_copy);
    }
  }
};

#endif