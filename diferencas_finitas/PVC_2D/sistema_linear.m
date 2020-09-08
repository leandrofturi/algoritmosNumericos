function [A,f] = sistema_linear(a,b,c,d,e,fun,n,m)

  N = n*m;
  A = sparse(N,N);
  f = zeros(N,1);

  % diagonal a_I
  for I = 1:N
    A(I,I) = a(I);
    f(I) = fun(I);
  endfor

  % diagonal b_I e c_I
  for I = 1:(N-1)
    A(I+1,I) = b(I+1);
    A(I,I+1) = c(I);
  endfor

  % diagonal d_I e e_I
  for I = 1:(N-n)
    A(I+n,I) = d(I+n);
    A(I,I+n) = e(I);
  endfor

endfunction