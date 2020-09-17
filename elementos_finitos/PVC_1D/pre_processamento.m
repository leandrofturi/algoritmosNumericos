function [n_nos, x, h, K, F, k_el, c_el, b_el, f_el] = pre_processamento(a, b, N_el, k_x, c_x, b_x, f_x);
  n_nos = N_el+1;
  h = 1/N_el;
  x(1, :) = linspace(a, b-h, n_nos-1);
  x(2, :) = linspace(a+h, b, n_nos-1);
  
  k_el(1, :) = k_x(x(1, :));
  k_el(2, :) = k_x(x(2, :));
  c_el(1, :) = c_x(x(1, :));
  c_el(2, :) = c_x(x(2, :));
  b_el(1, :) = b_x(x(1, :));
  b_el(2, :) = b_x(x(2, :));
  f_el(1, :) = f_x(x(1, :));
  f_el(2, :) = f_x(x(2, :));
  
  K = zeros(n_nos, n_nos);
  F = zeros(n_nos, 1);
endfunction
