clear
nome_matriz = "2-arc130.mat"

load(strcat("matrizes/", nome_matriz))
A = Problem.A;
printf("determinante de A = %d\n", det(A));
[n, n] = size(A);
printf("dimensao de A = %d %d\n", n, n);

b = A*ones(n, 1);

printf("\ndiagonal dominante: %s\n", diagonal_dominante(A));

w = [0:0.25:0.75, 1.25:0.25:2];
[MJ, MS, MSOR] = fatora(A, w);
[VJ, lambdaJ] = eig(MJ);
raioJ = max(abs(diag(lambdaJ)));
printf("\nraio espectral JACOBI: %.6f\n", raioJ);
[VS, lambdaS] = eig(MS);
raioS = max(abs(diag(lambdaS))); 
printf("raio espectral SEIDEL: %.6f\n", raioS);
for(i = 1:length(w))
  [VSOR, lambdaSOR] = eig(MSOR{i});
  raioSOR(i) = max(abs(diag(lambdaSOR)));
  printf("raio espectral SOR(%.2f): %.6f\n", w(i), raioSOR(i));
endfor;

tol = 10^-4;
nmaxiter = 500;
if(raioJ < 1)
  [xJ, erJ, iterJ] = jacobi(A, b, tol, nmaxiter);
endif;
if(raioS < 1)
  [xS, erS, iterS] = sor(A, b, tol, nmaxiter, 1);
endif;
for(i = 1:length(w))
  if(raioSOR(i) < 1)
    [xSOR{i}, erSOR{i}, iterSOR{i}] = sor(A, b, tol, nmaxiter, w(i));
  endif;
endfor;

[value, i] = min(raioSOR);
figure();
plot(1:iterJ, log(erJ), 1:iterS, log(erS), 1:iterSOR{i}, log(erSOR{i}))
xlabel("iteracao")
ylabel("log(erro)")
legend("Jacobi", "Seidel", sprintf("SOR(%.2f)", w(i)))
legend boxoff
% saveas(gcf, strcat("Q2-", nome_matriz, ".png"));