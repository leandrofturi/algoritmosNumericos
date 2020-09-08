function [u, tSOR, tDireto] = aplicacao2(tipo, w);
    
    if(tipo == "P")
      [u, tSOR, tDireto] = local_pvc2d(0,5000,0,1000,50,40,w);
    elseif(tipo == "G")
      [u, tSOR, tDireto] = local_pvc2d(0,5000,0,1000,500,400,w);
    else
      printf("Tamanho desconnhecido!");
    endif
  
endfunction

function [u, tSOR, tDireto] = local_pvc2d(a,b,c,d,n,m,w);

  hx = (b-a)/(n-1);
  hy = (d-c)/(m-1);
  x  = linspace(a,b,n);
  y  = linspace(c,d,m);

  [bx,by,gamma,f,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = local_aplicacao2(x,y,n,m);

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
  %local_grafico_velocidade(u,kappa,hx,hy,n,m)

endfunction

function [bx,by,gamma,f,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = local_aplicacao2(x,y,n,m)
  
  bx = repmat(0,1,n*m);
  by = repmat(0,1,n*m);
  gamma = repmat(0,1,n*m);
  f = repmat(0,1,n*m);
  %Localizacao dos pocos
  f(1, 600*m/1000 + n*(1500*n/5000 - 1)) = -250;
  f(1, 250*m/1000 + n*(3200*n/5000 - 1)) = -250;

  kappa = 1;

  left = 1;
  right = 1;
  bottom = 2;
  top = 2;
  
  gleft = repmat(100,n*m,1);
  gright = repmat(100,n*m,1);
  gbottom = repmat(0,n*m,1);
  gtop = repmat(0,n*m,1);
  
endfunction

function local_grafico_velocidade(u,kappa,hx,hy,n,m)
  
  figure();
  vx = zeros(m,n);
  vy = zeros(m,n);
  
  I = 1;
  vx(1,1) = (u(I+1) - u(I))/(hx);
  vy(1,1) = (u(I+n) - u(I))/(hy);
  I = (m-1)*n+1;
  vx(1,m) = (u(I+1) - u(I))/(hx);
  vy(1,m) = (u(I) - u(I-n))/(hy);
  I = n;
  vx(n,1) = (u(I) - u(I-1))/(hx);
  vy(n,1) = (u(I+n) - u(I))/(hy);
  I = n*m;
  vx(n,m) = (u(I) - u(I-1))/(hx);
  vy(n,m) = (u(I) - u(I-n))/(hy);
      
  % left
	for j = 2 : (m-1)
		for i = 1
      I = i + (j-1)*n;
      vx(j,i) = (u(I+1) - u(I))/(hx);
      vy(j,i) = (u(I+n) - u(I-n))/(2*hy);
		endfor
	endfor
  
  % right
	for j = 2 : (m-1)
		for i = n
      I = i + (j-1)*n;
			vx(j,i) = (u(I) - u(I-1))/(2);
      vy(j,i) = (u(I+n) - u(I-n))/(2*hy);
		endfor
	endfor

  % bottom
	for j = 1
		for i = 2 : (n-1)
      I = i + (j-1)*n;
			vx(j,i) = (u(I+1) - u(I-1))/(2*hx);
      vy(j,i) = (u(I+n) - u(I))/(2);
		endfor
	endfor

  % top
	for j = m
		for i = 2 : (n-1)
      I = i + (j-1)*n;
			vx(j,i) = (u(I+1) - u(I-1))/(2*hx);
      vy(j,i) = (u(I) - u(I-n))/(hy);
		endfor
	endfor

	for j = 2 : (m-1)
		for i = 2 : (n-1)
      I = i + (j-1)*n;
			vx(j,i) = (u(I+1) - u(I-1))/(2*hx);
      vy(j,i) = (u(I+n) - u(I-n))/(2*hy);
		endfor
	endfor
  
  vx = -kappa*vx;
  vy = -kappa*vy;
  
	quiver(vx,vy);
  xlabel('x');
  ylabel('y');
  
endfunction