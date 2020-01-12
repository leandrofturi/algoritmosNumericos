function [] = exe5();
  fileID = fopen('exe5.txt','w');

  %Q01
  %(a)
  fprintf(fileID, 'Q01(a)\n');
  a = 0;
  b = 10;
  ua = 40;
  ub = 200;
  for(n = [10, 50, 100])
    [p, q, r] = local_funcoesQ1(a, b, n);
    [x, u] = pvc(a, b, n, p, q, r, 1, ua, 0, 0, 0, 0, 1, ub, 0, 0, 0, 0);
    fprintf(fileID, 'n = %d\n', n);
    fprintf(fileID, "%f  ", u);
    fprintf(fileID, '\n');
    hf = figure();
    plot(x, u)
    xlabel ("x");
    ylabel ("T");
    legend("T(x)");
  endfor
  
  %(b)
  fprintf(fileID, '\nQ01(b)\n');
  a = 0;
  b = 10;
  ua = 40;
  sigma_b = 0;
  for(n = [10, 50, 100])
    [p, q, r] = local_funcoesQ1(a, b, n);
    [x, u] = pvc(a, b, n, p, q, r, 1, ua, 0, 0, 0, 0, 2, 0, sigma_b, 0, 0, 0);
    fprintf(fileID, 'n = %d\n', n);
    fprintf(fileID, "%f  ", u);
    fprintf(fileID, '\n');
    figure();
    plot(x, u)
    xlabel ("x");
    ylabel ("T");
    legend("T(x)");
  endfor

  %Q02
  fprintf(fileID, '\nQ02\n');
  colors = ['r', 'g', 'y', 'm', 'c', 'k'];
  for(n = [10, 50, 100])
    i = 1;
    figure();
    hold on;
    for(c_ref = [0.0001, 0.001, 0.01, 0.1])
      a = 0;
      b = 1;
      ua = 160;
      alfa_b = 0.001;
      beta_b = c_ref;
      gamma_b = c_ref*70;
      fprintf(fileID, 'c_ref = %d\n', c_ref);
      
      [p, q, r] = local_funcoesQ2(a, b, n, c_ref);
      [x, u] = pvc(a, b, n, p, q, r, 1, ua, 0, 0, 0, 0, 3, 0, 0, alfa_b, beta_b, gamma_b);
      fprintf(fileID, 'n = %d\n', n);
      fprintf(fileID, "%f  ", u);
      fprintf(fileID, '\n');
      plot(x, u, colors(i))
      i = i + 1;
    endfor
    xlabel ("x");
    ylabel ("u");
    legend("0.0001", "0.001", "0.01", "0.1");
    hold off;
  endfor
  
  fclose(fileID);
endfunction

function [p, q, r] = local_funcoesQ1(a, b, n);
  p = repmat(0, n, 1);
  q = repmat(-0.01, n, 1);
  r = repmat(-0.01*20, n, 1);
endfunction

function [p, q, r] = local_funcoesQ2(a, b, n, c_ref);
  x = linspace(a, b, n);
  
  p = repmat(0, n, 1);
  C = ((2*10 + 2*0.1)/(10*0.1))*c_ref;
  q = repmat(C/-0.001, n, 1);
  r = repmat(C*70/-0.001, n, 1);
endfunction

function [x, u] = pvc(a, b, n, p, q, r, tipo_a, ua, sigma_a, alfa_a, beta_a, gamma_a, tipo_b, ub, sigma_b, alfa_b, beta_b, gamma_b);
  h = (b-a)/(n-1);
  x = transpose(linspace(a, b, n));
  %[p, q, r] = funcoes(a, b, n);
  
  A = zeros(n, n);
  A(1, 1) = q(1) - (2/(h^2));
  A(1, 2) = (1/(h^2)) + (p(1)/(2*h));
  A(n, n-1) = (1/(h^2)) - (p(n)/(2*h));
  A(n, n) = q(n) - (2/(h^2));
  f = r;
  
  for i=2:(n-1)
    A(i, i-1) = (1/(h^2)) - (p(i)/(2*h)); %b_i
    A(i, i) = q(i) - (2/(h^2));           %a_i
    A(i, i+1) = (1/(h^2)) + (p(i)/(2*h)); %c_i
  endfor
  
  switch (tipo_a)
    case 1
      A(1, 1) = 1;
      A(1, 2) = 0;
      f(1) = ua;
    case 2
      A(1, 1) = A(1, 1) + (1/(h^2)) - (p(1)/(2*h));
      f(1) = f(1) + (((1/(h^2)) - (p(1)/(2*h)))*h*sigma_a);
    case 3
      A(1, 1) = A(1, 1) + (((1/(h^2)) - (p(1)/(2*h)))*(1 + (h*beta_a/alfa_a)));
      f(1) = f(1) + (((1/(h^2)) - (p(1)/(2*h)))*h*gamma_a/alfa_a);
  endswitch
  
  switch (tipo_b)
    case 1
      A(n, n-1) = 0;
      A(n, n) = 1;
      f(n) = ub;
    case 2
      A(n, n) = A(n, n) + (1/(h^2)) + (p(n)/(2*h));
      f(n) = f(n) - (((1/(h^2)) + (p(n)/(2*h)))*h*sigma_b);
    case 3
      A(n, n) = A(n, n) + (((1/(h^2)) + (p(n)/(2*h)))*(1 - (h*beta_b/alfa_b)));
      f(n) = f(n) - (((1/(h^2)) + (p(n)/(2*h)))*h*gamma_b/alfa_b);
  endswitch
  
  u = A\f;
  
endfunction