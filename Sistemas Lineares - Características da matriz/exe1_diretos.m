nome_matriz = "2-pivtol.mat"

load(strcat("matrizes/", nome_matriz))
A = Problem.A;
[n, n] = size(A);

[L, U, P] = lu(A);
printf("coef. nao nulos de A = %d\n", nnz(A));
figure();
spy(A);

printf("coef. nao nulos de L = %d\n", nnz(L));
figure();
spy(L);

printf("coef. nao nulos de U = %d\n", nnz(U));
figure();
spy(U);

printf("\n");

printf("taxa de preenchimento = %.6f\n", 100 - (nnz(A)/(nnz(L)+nnz(U)))*100);

printf("\n");

b = A*ones(n, 1);
x = A\b;

printf("distancia relativa entre a solucao exata e a aproximada = %.6f\n",
       norm(x, inf));

printf("\n");

printf("distancia relativa entre a matriz original e a matriz resultante da decomposicao LU = %.6f\n",
       norm(A-(P*L*U), inf)/norm(A, inf));

printf("\n");

printf("distancia relativa entre o vetor dos termos independentes original e o vetor resultante da decomposicao LU = %.6f\n",
       norm(b, inf));

printf("\n");

printf("norma do residuo = %.6f\n", norm(b - A*x, inf));

printf("\n");

printf("Numero de condicionamento da matriz = %.6f\n", cond(A));

printf("\n");