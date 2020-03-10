#include "tally_c.h"

// [[Rcpp::export]]
// Generate a list of frequencies
Rcpp::List calc_frequency_list(List rows, List cols, List col_names) {
  
  std::map<std::string, std::map<std::string, std::vector<double>>> row_map = convertToCMapping(rows);
  std::vector<std::vector<std::vector<double>>> col_data = convert_to_C_data(cols);
  std::vector<std::string> col_headers = get_headers(cols);
  std::map<std::string, std::vector<std::string>> inner_col_headers = get_inner_headers(cols);
  std::vector<std::string> df_col = convert_to_string_vector(col_names);
  std::vector<std::vector<std::vector<double>>> return_vector;

  List return_list;

  for(int i = 0; i < col_data.size(); i++) {
    return_list.push_back(apply_tally_c(col_data.at(i), inner_col_headers.at(col_headers.at(i)), df_col, i+1, row_map));
  }
  return return_list;
 
}


Rcpp::List apply_tally_c(std::vector<std::vector<double>> sub_col, 
                               std::vector<std::string> inner_col_headers,
                               std::vector<std::string> col_names,
                               int current_idx,
                               std::map<std::string, std::map<std::string, std::vector<double>>> row_map) {
  List return_list;

  std::vector<std::string> other_list_names = get_other_list_names(col_names, inner_col_headers);
  
  for(int i = 0; i < sub_col.size(); i++) {
    return_list.push_back(match_columns_c(sub_col.at(i), row_map.at(other_list_names.at(i)), current_idx));
  }
  
  
  return return_list;
}

Rcpp::List match_columns_c(std::vector<double> ref_vector, std::map<std::string, std::vector<double>> other_sublist, int current_idx) {
  
  Rcpp::List return_list;
  
  std::string current_idx_string = std::to_string(current_idx);
  std::vector<double> other_vector = other_sublist.at(current_idx_string);
  
  return_list.push_back(tabulate_c(ref_vector, other_vector));
  
  return return_list;
}

Rcpp::DataFrame tabulate_c(std::vector<double> ref, std::vector<double> other) {
  Rcpp::NumericVector return_vector;
  
  // Create a hashmap to count pairs
  std::map<std::pair<int, int>, int> pair_counts;
  
  for(int i = 0; i < ref.size(); i++) {
    std::pair<int, int> pair = std::make_pair(other.at(i), ref.at(i));
    
    if(pair_counts.count(pair) == 0) {
      // insert a new pair and init to 1 if it doesn't exist already
      pair_counts.insert(std::pair<std::pair<int,int>, int>(pair, 1));
    } else {
      // Increment existing pair
      pair_counts.at(pair) = pair_counts.at(pair) + 1;
    }
  }
  
  Rcpp::NumericVector other_vector;
  Rcpp::NumericVector ref_vector;
  Rcpp::NumericVector counts;
  
  // Count map should be finished, now generate a data frame
  for (auto& x: pair_counts) {
    other_vector.push_back(x.first.first);
    ref_vector.push_back(x.first.second);
    counts.push_back(x.second);
  }
  
  DataFrame df = DataFrame::create( Named("rvec") = other_vector, Named("cvec") = ref_vector, Named("Freq") = counts);
  return df;
}

// Helper functions below
Rcpp::NumericVector convert_to_numericalvector(std::vector<double> vec) {
  Rcpp::NumericVector return_vector;
  
  for(int i = 0; i < vec.size(); i++) {
    return_vector.push_back(vec.at(i));
  }
  
  return return_vector;
}

std::vector<std::string> convert_to_string_vector(List l) {
  std::vector<std::string> return_vector;
  for(int i =0; i < l.size(); i++) {
    return_vector.push_back(l.at(i));
  }
  
  return return_vector;
}

std::vector<std::string> get_other_list_names(std::vector<std::string> col_names, std::vector<std::string> inner_row_headers) {
  
  std::vector<std::string> return_vector;
  for(int i = 0; i< inner_row_headers.size(); i++) {
    return_vector.push_back(col_names.at(std::stoi(inner_row_headers.at(i))-1));
  }
  
  return return_vector;
}

std::vector<std::vector<std::vector<double>>> convert_to_C_data(List r_list) {
  
  std::vector<std::vector<std::vector<double>>> return_vector;
  std::vector<std::vector<double>> outer_vector;
  std::vector<double> inner_vector;
  
  for(int i = 0; i < r_list.size(); i++) {
    
    List r_outer_vector = Rcpp::as<List>(r_list[i]);
    
    for(int j = 0; j < r_outer_vector.size(); j++) {
      
      // Grab the inner vector of values
      inner_vector = Rcpp::as<std::vector<double>>(r_outer_vector[j]);
      
      // Insert inner vectors 
      outer_vector.push_back(inner_vector);
    }
    
    // Insert
    return_vector.push_back(outer_vector);
    outer_vector.clear();
  }
  return return_vector;
}

std::map<std::string, std::vector<std::string>> get_inner_headers(List r_list) {
  std::map<std::string, std::vector<std::string>> return_map;
  std::string key;
  std::vector<std::string> value;
  
  std::vector<std::string> outer_keys = Rcpp::as<std::vector<std::string>>(r_list.names());
  
  for(int i = 0; i < r_list.size(); i++) {
    List r_outer_vector = Rcpp::as<List>(r_list[i]);
    
    value = Rcpp::as<std::vector<std::string>>(r_outer_vector.names());
    key = outer_keys.at(i);
    
    return_map.insert(std::pair<std::string, std::vector<std::string>>(key,value));
  }
  return return_map;
}

std::vector<std::string> get_headers(List r_list) {
  std::vector<std::string> headers = Rcpp::as<std::vector<std::string>>(r_list.names());
  return headers;
}


std::map<std::string, std::map<std::string, std::vector<double>>> convertToCMapping(List r_list) { 
  
  // Map objects
  std::map<std::string, std::map<std::string, std::vector<double>>> return_map;
  std::map<std::string, std::vector<double>> inner_map;
  
  // inner vector
  std::vector<double> inner_vector;
  
  // Key names
  std::string outer_key;
  std::string inner_key;
  
  // Grab row/col name
  std::vector<std::string> outer_keys = Rcpp::as<std::vector<std::string>>(r_list.names());
  std::vector<std::string> inner_keys;
  
  for(int i = 0; i < r_list.size(); i++) {
    
    List r_outer_vector = Rcpp::as<List>(r_list[i]);
    
    // Get list of inner keys for that index
    inner_keys = Rcpp::as<std::vector<std::string>>(r_outer_vector.names());
    
    for(int j = 0; j < r_outer_vector.size(); j++) {
      
      // Grab the key for this index
      inner_key = inner_keys.at(j);
      
      // Grab the inner vector of values
      inner_vector = Rcpp::as<std::vector<double>>(r_outer_vector[j]);
      
      // Make a pair + insert
      inner_map.insert(std::make_pair(inner_key, inner_vector));
    }
    
    // Get outer key
    outer_key = outer_keys.at(i);
    
    // Make a pair + insert
    return_map.insert(std::make_pair(outer_key, inner_map));
    inner_map.clear();
  }
  
  return return_map;
}