function [x, iter, res] = gradientes_conjugados(A, b, tol, maxit)
  x = zeros(rows(A), 1);
  d = r = b;
  res(1, 1) = norm(r, inf);
  delta_new = delta_0 = norm(r, 2)^2;
  iter = 0;
  while(delta_new > tol^2*delta_0 && iter < maxit)
    iter = iter+1;
    v = A*d;
    lambda = delta_new/(transpose(d)*v);
    x = x + lambda*d;
    r = r - lambda*v;
    res(iter+1, 1) = norm(r, inf);
    delta_old = delta_new;
    delta_new = norm(r, 2)^2;
    beta = delta_new/delta_old;
    d = r + beta*d;
  endwhile
endfunction