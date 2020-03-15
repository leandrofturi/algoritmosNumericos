nome_matriz = "2-arc130.mat"

load(strcat("matrizes/", nome_matriz))
A = Problem.A;
printf("determinante de A = %d\n", det(A));
[n, n] = size(A);
printf("dimensao de A = %d %d\n", n, n);

b = A*ones(n, 1);

printf("\ndiagonal dominante: %s\n", diagonal_dominante(A));

w = [0:0.25:2];
[MJ, MS, MSOR] = fatora(A, w);
[VJ, lambdaJ] = eig(MJ);
printf("\nraio espectral JACOBI: %.6f\n", max(abs(diag(lambdaJ))));
[VS, lambdaS] = eig(MS);
printf("raio espectral SEIDEL: %.6f\n", max(abs(diag(lambdaS))));
for(i = 1:length(w))
  [VSOR, lambdaSOR] = eig(MSOR{i});
  printf("raio espectral SOR(%.2f): %.6f\n", w(i), max(abs(diag(lambdaSOR))));
endfor;

tol = 10^-6;
nmaxiter = 1000;
[x, er, iter] = jacobi(A, b, tol, nmaxiter);
for(i = w)
  [x, er, iter] = sor(A, b, tol, nmaxiter, i);
endfor;