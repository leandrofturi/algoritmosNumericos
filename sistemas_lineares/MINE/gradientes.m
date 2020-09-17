function [x, iter, res] = gradientes(A, b, tol, maxit)
  x = zeros(rows(A), 1);
  r = b;
  res(1, 1) = norm(r, inf);
  delta = delta_0 = norm(r, 2)^2;
  iter = 0;
  while(delta > tol^2*delta_0 && iter < maxit)
    iter = iter+1;
    v = A*r;
    lambda = delta/(transpose(r)*v);
    x = x + lambda*r;
    r = r - lambda*v;
    res(iter+1, 1) = norm(r, inf);
    delta = norm(r, 2)^2;
  endwhile
endfunction
