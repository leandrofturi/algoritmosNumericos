clear
nome_matriz = "4-poli_large.mat"

load(strcat("matrizes/", nome_matriz))
A = Problem.A;
printf("determinante de A = %d\n", det(A));
[n, n] = size(A);
printf("dimensao de A = %d %d\n", n, n);

warning('off', 'all');
[L, U, P] = lu(A);
printf("coef. nao nulos de A = %d\n", nnz(A));
figure();
spy(A);
% saveas(gcf, strcat("nnzA-", nome_matriz, ".png")); 
printf("coef. nao nulos de L = %d\n", nnz(L));
figure();
spy(L);
% saveas(gcf, strcat("nnzL-", nome_matriz, ".png"));
printf("coef. nao nulos de U = %d\n", nnz(U));
figure();
spy(U);
% saveas(gcf, strcat("nnzU-", nome_matriz, ".png"));

printf("taxa de preenchimento = %.6f\n", 100 - (nnz(A)/(nnz(L) + nnz(U)))*100);

b = A*ones(n, 1);
x = A\b;
dx = ones(n, 1) - x;
db = A*dx;
dA = A - (P*L*U);

printf("distancia relativa entre a solucao exata e a aproximada = %.6f\n",
       norm(dx, inf)/norm(x, inf));
printf("distancia relativa entre a matriz original e a matriz resultante da decomposicao LU = %.6f\n",
       norm(dA, inf)/norm(A, inf));
printf("distancia relativa entre o vetor dos termos independentes original e o vetor resultante da decomposicao LU = %.6f\n",
       norm(db, inf)/norm(b, inf));
printf("norma do residuo = %.6f\n", norm(b - A*x, inf));
printf("numero de condicionamento da matriz = %.6f\n", cond(A));