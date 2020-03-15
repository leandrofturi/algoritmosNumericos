nome_matriz = "2-arc130.mat"

load(strcat("matrizes/", nome_matriz))
A = Problem.A;
printf("determinante de A = %d\n", det(A));
[n, n] = size(A);
printf("dimensao de A = %d %d\n", n, n);

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

printf("\ntaxa de preenchimento = %.6f\n", 100 - (nnz(A)/(nnz(L)+nnz(U)))*100);

b = A*ones(n, 1);
x = A\b;

printf("\ndistancia relativa entre a solucao exata e a aproximada = %.6f\n",
       norm(x, inf)/norm(ones(n, 1), inf));

printf("\ndistancia relativa entre a matriz original e a matriz resultante da decomposicao LU = %.6f\n",
       norm(A-(P*L*U), inf)/norm(A, inf));

printf("\ndistancia relativa entre o vetor dos termos independentes original e o vetor resultante da decomposicao LU = %.6f\n",
       norm(b, inf));

printf("\nnorma do residuo = %.6f\n", norm(b - A*x, inf));

printf("\nnumero de condicionamento da matriz = %.6f\n", cond(A));