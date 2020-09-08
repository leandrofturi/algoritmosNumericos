function [n_nos, x, h, K, F, k_el, c_el, b_el, f_el] = pre_processamento(a, b, N_el, k_x, c_x, b_x, f_x);
  n_nos = N_el+1;
  x = linspace(a, b, n_nos);
  h = 1/N_el;
  
  k_el = k_x(x);
  c_el = c_x(x);
  b_el = b_x(x);
  f_el = f_x(x);
  
  K = zeros(n_nos, n_nos);
  F = zeros(n_nos, 1);
endfunction
