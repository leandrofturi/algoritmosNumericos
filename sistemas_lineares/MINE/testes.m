A = zeros(10, 10);
for(i = 1:10) A(i, i) = 10; endfor;
for(i = 2:10) A(i, i-1) = 1; endfor;
for(i = 1:9) A(i, i+1) = 1; endfor;
b = A*ones(10, 1);

tol = 10^-10;
maxit = 10^5;
[x, iter, res] = gradientes(A, b, tol, maxit);
[x, iter, res] = gradientes_conjugados(A, b, tol, maxit);
[x,flag,relres,iter,resvec] = pcg(A, b, tol, maxit);