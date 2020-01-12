function [x] = diagonal_dominante(A)
  [n, n] = size(A);
  x = 1;
  for i = 1:n
  aux = 0;
    for j = 1:(i-1)
      aux = aux + abs(A(i, j));
    endfor
    for j = (i+1):n
      aux = aux + abs(A(i, j));
    endfor
    if ((aux / abs(A(i, i))) >= 1)
      x = 0;
    endif
  endfor
endfunction