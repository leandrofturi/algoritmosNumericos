function [u, erro, tSOR, tDireto] = validacao2(tipo, w);
  
    if((tipo == "PP") && (size_equal(tipo, "PP")))
      [u, tSOR, tDireto] = local_pvc2d(0,1,0,1,5,5,w);
      m = n = 5;
    elseif((tipo == "P") && (size_equal(tipo, "P")))
      [u, tSOR, tDireto] = local_pvc2d(0,1,0,1,50,50,w);
      m = n = 50;
    elseif((tipo == "G") && (size_equal(tipo, "G")))
      [u, tSOR, tDireto] = local_pvc2d(0,1,0,1,100,250,w);
      m = 250; n = 100;
    elseif((tipo == "GG") && (size_equal(tipo, "GG")))
      [u, tSOR, tDireto] = local_pvc2d(0,1,0,1,500,500,w);
      m = n = 500;
    else
      printf("Tamanho desconnhecido!");
    endif
  
  x  = linspace(0,1,n);
  y  = linspace(0,1,m);
  U = zeros(n*m,1);
  I = 1;
  for j = 1:m
  	for i = 1:n
      U(I,1) = 10*x(i)*y(j)*(1-x(i))*(1-y(j))*exp(x(i)^4.5);
      I ++;
    endfor
  endfor
  %grafico_solucao(U,x,y,n,m)
  erro = norm(U-u,inf);
  
endfunction

function [u, tSOR, tDireto] = local_pvc2d(a,b,c,d,n,m,w);

  hx = (b-a)/(n-1);
  hy = (d-c)/(m-1);
  x  = linspace(a,b,n);
  y  = linspace(c,d,m);

  [bx,by,gamma,f,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = local_validacao2(x,y,n,m);

  [ai,bi,ci,di,ei] = coeficientes(hx,hy,kappa,bx,by,gamma,n,m);

  [ai,bi,ci,di,ei,f] = condicoes_contorno(ai,bi,ci,di,ei,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,hx,hy,kappa);
  
  [A,f] = sistema_linear(ai,bi,ci,di,ei,f,n,m);
  
  tDireto = tic();
  u = A\f;
  tDireto = toc(tDireto);
  
  tSOR = tic();
  [u,iter,er] = sor(ai,bi,ci,di,ei,f,n,m,w,0.00001,1000);
  tSOR = toc(tSOR);

  %grafico_solucao(u,x,y,n,m)

endfunction

function [bx,by,gamma,f,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = local_validacao2(x,y,n,m)

  I = 1;
  for j = 1:m
  	for i = 1:n 
      bx(I) = 1.0;
      by(I) = 20*y(j);
      gamma(I) = 1.0;
      f(I) = -((337.5*(exp(x(i)^4.5))*(x(i)^4.5 - 0.733333*x(i)^3.5 + 0.6*(x(i)^9) - 0.6*(x(i)^8) + 0.059259)*(y(j)-1)*y(j)) + (20*(exp(x(i)^4.5))*(x(i)-1)*x(i)) + bx(I)*(-45*(exp(x(i)^4.5))*(-x(i)^5.5 + x(i)^4.5 - 0.444444*x(i) + 0.222222)*(y(j)-1)*y(j)) + by(I)*(10*(exp(x(i)^4.5))*(x(i)-1)*x(i)*(2*y(j)-1)) + gamma(I)*(10*x(i)*y(j)*(1-x(i))*(1-y(j))*exp(x(i)^4.5)));
      I++;
    endfor
  endfor

  kappa = 1.0;

  left = 1;
  right = 1;
  bottom = 1;
  top = 1;
  
  gleft = repmat(0,n*m,1);
  gright = repmat(0,n*m,1);
  gbottom = repmat(0,n*m,1);
  gtop = repmat(0,n*m,1);
  
  % Solucao real e suas derivadas
  
  %10*x*y*(1-x)*(1-y)*e^(x^4.5)
  
  %-45*(e^(x^4.5))*(-x^5.5 + x^4.5 - 0.444444*x + 0.222222)*(y-1)*y
  %337.5*(e^(x^4.5))*(x^4.5 - 0.733333*x^3.5 + 0.6*(x^9) - 0.6*(x^8) + 0.059259)*(y-1)*y
  
  %10*(e^(x^4.5))*(x-1)*x*(2*y-1)
  %20*(e^(x^4.5))*(x-1)*x
  
endfunction
