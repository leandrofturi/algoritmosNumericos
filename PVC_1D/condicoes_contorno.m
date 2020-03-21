function [x, u] = condicoes_contorno(a, b, n, p, q, r, tipo_a, ua, sigma_a, alfa_a, beta_a, gamma_a, tipo_b, ub, sigma_b, alfa_b, beta_b, gamma_b);

  h = (b-a)/(n-1);
  x = transpose(linspace(a, b, n));
  
  A = zeros(n, n);
  A(1, 1) = q(1) - (2/(h^2));
  A(1, 2) = (1/(h^2)) + (p(1)/(2*h));
  A(n, n-1) = (1/(h^2)) - (p(n)/(2*h));
  A(n, n) = q(n) - (2/(h^2));
  f = r;
  
  for i=2:(n-1)
    A(i, i-1) = (1/(h^2)) - (p(i)/(2*h)); %b_i
    A(i, i) = q(i) - (2/(h^2));           %a_i
    A(i, i+1) = (1/(h^2)) + (p(i)/(2*h)); %c_i
  endfor
  
  switch (tipo_a)
    case 1
      A(1, 1) = 1;
      A(1, 2) = 0;
      f(1) = ua;
    case 2
      A(1, 1) = A(1, 1) + (1/(h^2)) - (p(1)/(2*h));
      f(1) = f(1) + (((1/(h^2)) - (p(1)/(2*h)))*h*sigma_a);
    case 3
      A(1, 1) = A(1, 1) + (((1/(h^2)) - (p(1)/(2*h)))*(1 + (h*beta_a/alfa_a)));
      f(1) = f(1) + (((1/(h^2)) - (p(1)/(2*h)))*h*gamma_a/alfa_a);
  endswitch
  
  switch (tipo_b)
    case 1
      A(n, n-1) = 0;
      A(n, n) = 1;
      f(n) = ub;
    case 2
      A(n, n) = A(n, n) + (1/(h^2)) + (p(n)/(2*h));
      f(n) = f(n) - (((1/(h^2)) + (p(n)/(2*h)))*h*sigma_b);
    case 3
      A(n, n) = A(n, n) + (((1/(h^2)) + (p(n)/(2*h)))*(1 - (h*beta_b/alfa_b)));
      f(n) = f(n) - (((1/(h^2)) + (p(n)/(2*h)))*h*gamma_b/alfa_b);
  endswitch
  
  u = A\f;
  
endfunction