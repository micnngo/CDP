import numpy as np
import dt2b

def fit_dt2b(filename, Kr, Kc, n = 10, alpha = 0.5, beta_row = 0.5, beta_column = 0.5, verbose = True):
  data = np.genfromtxt(filename, delimiter=',', dtype=bytes).astype(str)
  col_w_z,row_a_z,joint,idx2vals,vals2idx = dt2b.dt2b(data, Kr, Kc, n, alpha, beta_row, beta_column, it = 500, verbose = verbose)
  return({'joint': joint, 'phi_r': row_a_z, 'phi_c': col_w_z, 'idx2vals': idx2vals})


