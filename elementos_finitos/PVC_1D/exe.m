a = 0;
b = 1;
N_el = 4;

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
gama_a = 1;
alfa_b = 0;
beta_b = 1;
gama_b = 1;

elem = 2.^(3:8);
j = 1;
for i = elem
  [n_nos, x, h, K, F, k_el, c_el, b_el, f_el] = pre_processamento(a, b, i, @k_x, @c_x, @b_x, @f_x);
  [K, F, u] = processamento(i, n_nos, x, K, F, alfa_a, beta_a, gama_a, alfa_b, beta_b, gama_b, k_el, c_el, b_el, f_el);
  x = linspace(a, b, n_nos);
  s = transpose(x.^3 - x + 1);
  er(j) = pos_processamento(x, s, u);
  j = j+1;
endfor

x = log(1./elem);
y = log(er);
p = polyfit(x, y, 1)
plot(x, y, ".", x, p(1)*x+p(2))