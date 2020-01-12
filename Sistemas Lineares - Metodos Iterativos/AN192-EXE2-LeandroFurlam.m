function [xJ, xS, xSOR] = AN192-EXE4-LeandroFurlam(arq, w, tol, nmaxiter, nome)
  source jacobi.m;
  source sor.m;
  source fatora.m;
  source diagonal_dominante.m;

  % Abrindo o arquivo
  fileID = fopen(nome,'w');
  fprintf(fileID,'%s\n', arq);
  
  % Carregando a matriz
  load (arq);
  A = Problem.A;
  [n, n] = size(A);
  fprintf(fileID,'Dimensao: %d\n', n);
  
  % Vetor dos termos independentes
  b = A*ones(n, 1);
  
  % Verificar se e diagonal dominante
  ehDiagonal = diagonal_dominante(A);
  fprintf(fileID,'Diagonal dominante: %d\n', ehDiagonal);
  
  % Calculo do raio espectral
  [MJ, MS, MSOR] = fatora(A, w);
  
  lambdaJ = eig(MJ);
  for i = 1:rows(lambdaJ)
    lambdaJ(i) = norm(lambdaJ(i), 2);
  endfor
  raioMJ = max(abs(lambdaJ));
  fprintf(fileID,'\nRaio Espectral\n');
  fprintf(fileID,'Jacobi: %d\n', raioMJ);
  
  lambdaS = eig(MS);
  for i = 1:rows(lambdaS)
    lambdaS(i) = norm(lambdaS(i), 2);
  endfor
  raioMS = max(abs(lambdaS));
  fprintf(fileID,'Seidel: %d\n', raioMS);
  
  lambdaSOR = eig(MSOR);
  for i = 1:rows(lambdaSOR)
    lambdaSOR(i) = norm(lambdaSOR(i), 2);
  endfor
  raioMSOR = max(abs(lambdaSOR));
  fprintf(fileID,'SOR: %d\n', raioMSOR);
  
  % Calculo dos metodos
  fprintf(fileID,'\nMetodos\n');
  if (raioMJ < 1)
    [xJ, erJ, iterJ] = jacobi(A, b, tol, nmaxiter);
  else
    xJ = 0;
    erJ = 0;
    iterJ = 0;
  endif
  fprintf(fileID,'Jacobi: %d\n', iterJ);
  fprintf(fileID,'%f\n', erJ);
  fprintf(fileID,'\n');
  
  if (raioMS < 1)
    [xS, erS, iterS] = sor(A, b, tol, nmaxiter, 1);
  else
    xS = 0;
    erS = 0;
    iterS = 0;
  endif
  fprintf(fileID,'Seidel: %d\n', iterS);
  fprintf(fileID,'%f\n', erS);
  fprintf(fileID,'\n');
  
  if (raioMSOR < 1)
    [xSOR, erSOR, iterSOR] = sor(A, b, tol, nmaxiter, w);
  else
    xSOR = 0;
    erSOR = 0;
    iterSOR = 0;
  endif
  fprintf(fileID,'SOR: %d\n', iterSOR);
  fprintf(fileID,'%f\n', erSOR);
  fprintf(fileID,'\n');
  
  % Fechando o arquivo
  fclose(fileID);
  
  % raioMJ, raioMS, raioMSOR, tol, nmaxiter, xJ, erJ, iterJ, xS, erS, iterS, xSOR, erSOR, iterSOR
endfunction