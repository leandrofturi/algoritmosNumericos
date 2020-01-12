function [x, er, iter] = sor(A, b, tol, nmaxiter, w)
  [n, n] = size(A);
  x = x_seidel = zeros(n, 1);
  iter = 1;
  er(iter) = Inf;
  
  while ((iter <= nmaxiter) && (tol <= er(iter)))
    x0 = x;
    x0_seidel = x_seidel;
    for i = 1:n
      aux = 0;
      for j = 1:(i-1)
         aux = aux + A(i, j)*x_seidel(j);
      endfor
      for j = (i+1):n
         aux = aux + A(i, j)*x0_seidel(j);
      endfor
      x_seidel(i) = (1/A(i, i))*(b(i) - aux);
      endfor
    x = w*x_seidel + (1-w)*x0;
    
    iter = iter + 1;
    er(iter) = norm(x - x0, inf) / norm(x, inf);
  endwhile
endfunction