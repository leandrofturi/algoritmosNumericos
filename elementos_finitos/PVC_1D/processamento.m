function [K, F, u] = processamento(n_nos, h, x, alfa_a, beta_a, gamma_a, alfa_b, beta_b, gamma_b, k_el, c_el, b_el, f_el, K, F);
  for i = 1:(n_nos-1)
    K(i, i) = K(i, i) + k_int(k_el, x, i) + c_int(c_el, x, i, x(i+1)) + b_int(b_el, x, i, x(i+1), x(i+1));
    K(i, i+1) = K(i, i+1) - k_int(k_el, x, i) - c_int(c_el, x, i, x(i+1)) - b_int(b_el, x, i, x(i+1), x(i));
    K(i+1, i) = K(i+1, i) - k_int(k_el, x, i) - c_int(c_el, x, i, x(i)) - b_int(b_el, x, i, x(i), x(i+1));
    K(i+1, i+1) = K(i+1, i+1) + k_int(k_el, x, i) + c_int(c_el, x, i, x(i)) + b_int(b_el, x, i, x(i), x(i));
    
    F(i) = F(i) - f_int(f_el, x, i, x(i+1));
    F(i+1) = F(i+1) + f_int(f_el, x, i, x(i));
  endfor
  
  % Condição de contorno de Dirichlet em x = a
  if(alfa_a == 0)
    K(1, :) = 0;
    K(1, 1) = 1;
    F(1) = gamma_a/beta_a;

  % Condição de contorno de Newmann ou mista em x = a
  else
    K(1, 1) = K(1, 1) - k_el(1)*beta_a/alfa_a;
    F(1) = F(1) - k_el(1)*gamma_a/alfa_a;
  endif

  % Condição de contorno de Dirichlet em x = b
  if(alfa_b == 0)
    K(n_nos, :) = 0;
    K(n_nos, n_nos) = 1;
    F(n_nos) = gamma_b/beta_b;

  % Condição de contorno de Newmann ou mista em x = b
  else
    K(n_nos, n_nos) = K(n_nos, n_nos) + k_el(n_nos)*beta_b/alfa_b;
    F(n_nos) = F(n_nos) + k_el(n_nos)*gamma_b/alfa_b;
  endif
  
  u = K\F;
endfunction

function y = k_int(k_el, x, i)
  h = x(i+1) - x(i);
  x1 = x(i);
  x2 = x(i+1);
  k1 = k_el(i);
  k2 = k_el(i+1);
  y = (1/h^3)*( ( (x2^2/2)*(-k1 + k2) + x2*(k1*x2 - k2*x1) ) -
                ( (x1^2/2)*(-k1 + k2) + x1*(k1*x2 - k2*x1) ) );
endfunction

function y = c_int(c_el, x, i, A)
  h = x(i+1) - x(i);
  x1 = x(i);
  x2 = x(i+1);
  c1 = c_el(i);
  c2 = c_el(i+1);
  y = (1/h^3)*( ( (x2^3/3 - (x2^2/2)*A)*(-c1 + c2) + (x2^2/2 - A*x2)*(c1*x2 - c2*x1) ) - 
                ( (x1^3/3 - (x1^2/2)*A)*(-c1 + c2) + (x1^2/2 - A*x1)*(c1*x2 - c2*x1) ) );
endfunction

function y = b_int(b_el, x, i, A, B)
  h = x(i+1) - x(i);
  x1 = x(i);
  x2 = x(i+1);
  b1 = b_el(i);
  b2 = b_el(i+1);
  y = (1/h^3)*( ( (x2^4/4 - (x2^3/3)*(A+B) + (x2^2/2)*A*B)*(-b1 + b2) + (x2^3/3 - (x2^2/2)*(A+B) + x2*A*B)*(b1*x2 - b2*x1) ) - 
                ( (x1^4/4 - (x1^3/3)*(A+B) + (x1^2/2)*A*B)*(-b1 + b2) + (x1^3/3 - (x1^2/2)*(A+B) + x1*A*B)*(b1*x2 - b2*x1) ) );
endfunction

function y = f_int(f_el, x, i, A)
  h = x(i+1) - x(i);
  x1 = x(i);
  x2 = x(i+1);
  f1 = f_el(i);
  f2 = f_el(i+1);
  y = (1/h^2)*( ( (x2^3/3 - (x2^2/2)*A)*(-f1 + f2) + (x2^2/2 - A*x2)*(f1*x2 - f2*x1) ) - 
                ( (x1^3/3 - (x1^2/2)*A)*(-f1 + f2) + (x1^2/2 - A*x1)*(f1*x2 - f2*x1) ) );
endfunction
