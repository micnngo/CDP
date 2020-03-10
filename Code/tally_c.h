#ifndef TALLY_C_H
#define TALLY_C_H

#include <Rcpp.h>
#include <RcppParallel.h>

#include <map>
#include <vector>
#include <string>
using namespace Rcpp;

Rcpp::List calc_frequency_list(List rows, List cols, List col_names);
std::map<std::string, std::map<std::string, std::vector<double>>> convertToCMapping(List r_list);
std::vector<std::vector<std::vector<double>>> convert_to_C_data(List r_list);
std::vector<std::string> get_headers(List r_list);
Rcpp::List apply_tally_c(std::vector<std::vector<double>> sub_col, 
                               std::vector<std::string> inner_col_headers,
                               std::vector<std::string> col_names,
                               int current_idx,
                               std::map<std::string, std::map<std::string, std::vector<double>>> row_map);

std::vector<std::string> convert_to_string_vector(List l);
std::map<std::string, std::vector<std::string>> get_inner_headers(List r_list);
std::vector<std::string> get_other_list_names(std::vector<std::string> col_names, std::vector<std::string> inner_row_headers);
Rcpp::List match_columns_c(std::vector<double> ref_vector, std::map<std::string, std::vector<double>> other_sublist, int current_idx);
Rcpp::NumericVector convert_to_numericalvector(std::vector<double> vec);
Rcpp::DataFrame tabulate_c(std::vector<double> ref, std::vector<double> other);

#endif