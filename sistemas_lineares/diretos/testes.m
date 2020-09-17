clear
[fname, fpath, fltidx] = uigetfile("*.mat", "Escolha a matriz");
nome_matriz = strrep(fname, ".mat", "");

diary(strrep(strcat(fpath, fname, "-diretos"), ".mat", ".txt"))

printf("%s\n", nome_matriz);
load(strcat(fpath, fname));
A = Problem.A;
printf("determinante de A = %.6f\n", det(A));
[n, n] = size(A);
printf("dimensao de A = %d %d\n", n, n);
printf("\n");

[L, U, P] = lu(A);
printf("coef. nao nulos de A = %d\n", nnz(A));
hold on;
spy(A);
printf("coef. nao nulos de L = %d\n", nnz(L));
figure();
spy(L);
printf("coef. nao nulos de U = %d\n", nnz(U));
figure();
spy(U);

printf("taxa de preenchimento = %.4f\n", 100 - (nnz(A)/(nnz(L) + nnz(U)))*100);
printf("\n");

b = A*ones(n, 1);
x = A\b;
dx = ones(n, 1) - x;
db = A*dx;
dA = A - (P*L*U);

printf("distancia relativa entre a solucao exata e a aproximada = %.4f\n",
       norm(dx, inf)/norm(x, inf));
printf("distancia relativa entre a matriz original e a matriz resultante da decomposicao LU = %.4f\n",
       norm(dA, inf)/norm(A, inf));
printf("distancia relativa entre o vetor dos termos independentes original e o vetor resultante da decomposicao LU = %.4f\n",
       norm(db, inf)/norm(b, inf));
printf("norma do residuo = %.4f\n", norm(b - A*x, inf));
printf("numero de condicionamento da matriz = %.4f\n", cond(A));

diary off