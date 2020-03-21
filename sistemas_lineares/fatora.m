function[MJ, MS, MSOR] = fatora(A, w)

  [n, n] = size(A);
  D = diag(diag(A));
  E = tril(A, -1);
  F = triu(A, 1);
  
  MJ = (-1)*inv(D)*(E + F);
  MS = (-1)*inv(E + D)*F;
  for(i = 1:length(w))
    MSOR{i} = inv(D + w(i)*E)*((1-w(i))*D - w(i)*F);
  endfor;

endfunction