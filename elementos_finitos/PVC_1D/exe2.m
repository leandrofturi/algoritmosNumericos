a = 0;
b = 1;
N_el = 5;

function y = k_x(x);
  y = sin(x);
endfunction

function y = c_x(x);
  y = x;
endfunction

function y = b_x(x);
  y = cos(x);
endfunction

function y = f_x(x);
  for i=1:size(x)(:,2)
    y(i) = x(i)^2*(cos(x(i))+2)-2*x(i)*cos(x(i))+cos(x(i))-2*sin(x(i));
  endfor
endfunction

alfa_a = 0;
beta_a = 1;
gama_a = 1;
alfa_b = 1;
beta_b = -1;
gama_b = 0;

[n_nos, x, h, K, F, k_el, c_el, b_el, f_el] = pre_processamento(a, b, N_el, @k_x, @c_x, @b_x, @f_x);
[K, F, u] = processamento(N_el, n_nos, x, K, F, alfa_a, beta_a, gama_a, alfa_b, beta_b, gama_b, k_el, c_el, b_el, f_el);

x = linspace(a, b, n_nos);
s = transpose(x.^2 + 1);
pos_processamento(x, s, u);