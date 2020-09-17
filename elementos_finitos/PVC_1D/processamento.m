function [K, F, u] = processamento(N_el, n_nos, x, K, F, alfa_a, beta_a, gama_a, alfa_b, beta_b, gama_b, k_el, c_el, b_el, f_el);
  for i = 1:N_el
    K(i:(i+1), i:(i+1)) = K(i:(i+1), i:(i+1)) + K_el(x(:, i), k_el(:, i), c_el(:, i), b_el(:, i));
    F(i:(i+1)) = F(i:(i+1)) + F_el(x(:, i), f_el(:, i));
  endfor
  
  % Condição de contorno de Dirichlet em x = a
  if(alfa_a == 0)
    F(1) = gama_a/beta_a;
    F(2) = F(2) - K(2, 1)*gama_a/beta_a;
    K(1, :) = K(:, 1) = 0;
    K(1, 1) = 1;
  % Condição de contorno de Newmann ou mista em x = a
  else
    K(1, 1) = K(1, 1) - k_el(1)*beta_a/alfa_a;
    F(1) = F(1) - k_el(1)*gama_a/alfa_a;
  endif

  % Condição de contorno de Dirichlet em x = b
  if(alfa_b == 0)
    F(n_nos) = gama_b/beta_b;
    F(n_nos-1) = F(n_nos-1) - K(n_nos-1, n_nos)*gama_b/beta_b;
    K(n_nos, :) = K(:, n_nos) = 0;
    K(n_nos, n_nos) = 1;
  % Condição de contorno de Newmann ou mista em x = b
  else
    K(n_nos, n_nos) = K(n_nos, n_nos) + k_el(n_nos)*beta_b/alfa_b;
    F(n_nos) = F(n_nos) + k_el(n_nos)*gama_b/alfa_b;
  endif
  
  u = K\F;
endfunction

function K = K_el(x_el, k_el, c_el, b_el)
  K(1, 1) = k_int(k_el, x_el) + c_int(c_el, x_el, x_el(2)) + b_int(b_el, x_el, x_el(2), x_el(2));
  K(1, 2) = - k_int(k_el, x_el) - c_int(c_el, x_el, x_el(2)) - b_int(b_el, x_el, x_el(2), x_el(1));
  K(2, 1) = - k_int(k_el, x_el) - c_int(c_el, x_el, x_el(1)) - b_int(b_el, x_el, x_el(1), x_el(2));
  K(2, 2) = k_int(k_el, x_el) + c_int(c_el, x_el, x_el(1)) + b_int(b_el, x_el, x_el(1), x_el(1));
endfunction

function F = F_el(x_el, f_el)
  F(1, 1) = - f_int(f_el, x_el, x_el(2));
  F(2, 1) = f_int(f_el, x_el, x_el(1));
endfunction

function y = k_int(k_el, x_el)
  h = x_el(2) - x_el(1);
  x1 = x_el(1);
  x2 = x_el(2);
  k1 = k_el(2);
  k2 = k_el(2);
  y = (1/h^3)*( ( (x2^2/2)*(-k1 + k2) + x2*(k1*x2 - k2*x1) ) -
                ( (x1^2/2)*(-k1 + k2) + x1*(k1*x2 - k2*x1) ) );
endfunction

function y = c_int(c_el, x_el, A)
  h = x_el(2) - x_el(1);
  x1 = x_el(1);
  x2 = x_el(2);
  c1 = c_el(1);
  c2 = c_el(2);
  y = (1/h^3)*( ( (x2^3/3 - (x2^2/2)*A)*(-c1 + c2) + (x2^2/2 - A*x2)*(c1*x2 - c2*x1) ) - 
                ( (x1^3/3 - (x1^2/2)*A)*(-c1 + c2) + (x1^2/2 - A*x1)*(c1*x2 - c2*x1) ) );
endfunction

function y = b_int(b_el, x_el, A, B)
  h = x_el(2) - x_el(1);
  x1 = x_el(1);
  x2 = x_el(2);
  b1 = b_el(1);
  b2 = b_el(2);
  y = (1/h^3)*( ( (x2^4/4 - (x2^3/3)*(A+B) + (x2^2/2)*A*B)*(-b1 + b2) + (x2^3/3 - (x2^2/2)*(A+B) + x2*A*B)*(b1*x2 - b2*x1) ) - 
                ( (x1^4/4 - (x1^3/3)*(A+B) + (x1^2/2)*A*B)*(-b1 + b2) + (x1^3/3 - (x1^2/2)*(A+B) + x1*A*B)*(b1*x2 - b2*x1) ) );
endfunction

function y = f_int(f_el, x_el, A)
  h = x_el(2) - x_el(1);
  x1 = x_el(1);
  x2 = x_el(2);
  f1 = f_el(1);
  f2 = f_el(2);
  y = (1/h^2)*( ( (x2^3/3 - (x2^2/2)*A)*(-f1 + f2) + (x2^2/2 - A*x2)*(f1*x2 - f2*x1) ) - 
                ( (x1^3/3 - (x1^2/2)*A)*(-f1 + f2) + (x1^2/2 - A*x1)*(f1*x2 - f2*x1) ) );
endfunction
