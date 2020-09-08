a = 0;
b = 10;
ua = 40;
ub = 200;
for(n = [10, 50, 100])

  p = repmat(0, n, 1);
  q = repmat(-0.01, n, 1);
  r = repmat(-0.01*20, n, 1);

  [x, u] = condicoes_contorno(a, b, n, p, q, r, 1, ua, 0, 0, 0, 0, 1, ub, 0, 0, 0, 0);
  figure();
  plot(x, u)
  xlabel("x");
  ylabel("T");
  legend("T(x)");
endfor


a = 0;
b = 10;
ua = 40;
sigma_b = 0;
for(n = [10, 50, 100])

  p = repmat(0, n, 1);
  q = repmat(-0.01, n, 1);
  r = repmat(-0.01*20, n, 1);

  [x, u] = condicoes_contorno(a, b, n, p, q, r, 1, ua, 0, 0, 0, 0, 2, 0, sigma_b, 0, 0, 0);
  figure();
  plot(x, u)
  xlabel("x");
  ylabel("T");
  legend("T(x)");
endfor


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
     
    x = linspace(a, b, n);
  
    p = repmat(0, n, 1);
    C = ((2*10 + 2*0.1)/(10*0.1))*c_ref;
    q = repmat(C/-0.001, n, 1);
    r = repmat(C*70/-0.001, n, 1);

    [x, u] = condicoes_contorno(a, b, n, p, q, r, 1, ua, 0, 0, 0, 0, 3, 0, 0, alfa_b, beta_b, gamma_b);
    plot(x, u, colors(i))
    i = i + 1;
  endfor
  xlabel ("x");
  ylabel ("u");
  legend("0.0001", "0.001", "0.01", "0.1");
  hold off;
endfor