function [x,iter,er] = sor(a,b,c,d,e,f,n,m,w,tol,maxiter);
 
  x0 = x = zeros(n*m, 1);
  iter = 1;
  er(iter) = Inf;
   
  A = zeros(n*m, 1);
  B = zeros(n*m, 1);
  C = zeros(n*m, 1);
  D = zeros(n*m, 1);
  E = zeros(n*m, 1);
  F = zeros(n*m, 1);
   
  while ((iter <= maxiter) && (tol <= er(iter)))
    x0 = x;
     
    I = 1;
      A(I) = a(I);
      C(I) = c(I)*x0(I+1);
      E(I) = e(I)*x0(I+n);
      F(I) = f(I);
    for I = 2:n
      B(I) = b(I)*(x(I-1));
      A(I) = a(I);
      C(I) = c(I)*x0(I+1);
      E(I) = e(I)*x0(I+n);
      F(I) = f(I);
    endfor
    for I = (n+1):((m-1)*n)
      D(I) = d(I)*x0(I-n);
      B(I) = b(I)*(x(I-1));
      A(I) = a(I);
      C(I) = c(I)*x0(I+1);
      E(I) = e(I)*x0(I+n);
      F(I) = f(I);
    endfor
    for I = ((m-1)*n+1):(m*n-1)
      D(I) = d(I)*x0(I-n);
      B(I) = b(I)*(x(I-1));
      A(I) = a(I);
      C(I) = c(I)*x0(I+1);
      F(I) = f(I);
    endfor
    I = n*m;
      D(I) = d(I)*x0(I-n);
      B(I) = b(I)*(x(I-1));
      A(I) = a(I);
      F(I) = f(I);
     
    x = w./A.*(F - D - B - C - E) + (1-w).*x0;
     
    iter++;
    er(iter) = norm(x-x0,inf)/norm(x,inf);
  endwhile
  iter--;
 
endfunction