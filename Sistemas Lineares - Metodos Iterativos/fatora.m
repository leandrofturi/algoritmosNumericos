function [MJ, MS, MSOR] = fatora(A, w)
  D = diag(diag(A));
  E = tril(A, -1);
  F = tril(A, 1);
  
  MJ = (-1)*inv(D)*(E + F);
  MS = (-1)*inv(E + D)*F;
  MSOR = inv(D + w*E)*((1-w)*D - w*F);
endfunction