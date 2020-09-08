a = 0;
b = 1;
N_el = 5;

function y = k_x(x);
  y = x;
endfunction

function y = c_x(x);
  y = ones(1, length(x));
endfunction

function y = b_x(x);
  y = x.^2;
endfunction

function y = f_x(x);
  y = x.^5 - x.^3 - 5*x.^2;
endfunction

alfa_a = 0;
beta_a = 1;
gamma_a = 1;
alfa_b = 2;
beta_b = 1;
gamma_b = 5;

[n_nos, x, h, K, F, k_el, c_el, b_el, f_el] = pre_processamento(a, b, N_el, @k_x, @c_x, @b_x, @f_x);
[K, F, u] = processamento(n_nos, h, x, alfa_a, beta_a, gamma_a, alfa_b, beta_b, gamma_b, k_el, c_el, b_el, f_el, K, F);

s = transpose(x.^3 - x + 1);
pos_processamento(x, s, u);