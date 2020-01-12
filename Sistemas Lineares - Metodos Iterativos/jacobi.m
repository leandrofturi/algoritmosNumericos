function [x, er, iter] = jacobi(A, b, tol, nmaxiter)
  [n, n] = size(A);
  x = zeros(n, 1);
  iter = 1;
  er(iter) = Inf;
  
  while ((iter <= nmaxiter) && (tol <= er(iter)))
    x0 = x;
    for i = 1:n
      aux = 0;
      for j = 1:(i-1)
         aux = aux + A(i, j)*x0(j);
      endfor
      for j = (i+1):n
         aux = aux + A(i, j)*x0(j);
      endfor
      x(i) = (1/A(i, i))*(b(i) - aux);      
    endfor
    
    iter = iter + 1;
    er(iter) = norm(x - x0, inf) / norm(x, inf);
  endwhile
endfunction

%A = [6, -2, 0, -1; 1, 5, -2, 0; 0, 2, 6, -2; 0, 0, -1, 4]
%b = [3;6;4;3]