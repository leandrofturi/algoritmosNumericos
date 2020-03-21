function [u, tSOR, tDireto] = validacao1(tipo, w);
  
    if(tipo == "PP")
      [u, tSOR, tDireto] = local_pvc2d(0,1,0,1,5,6,w);
    else
      printf("Tamanho desconnhecido!");
    endif
  
endfunction

function [u, tSOR, tDireto] = local_pvc2d(a,b,c,d,n,m,w);

  hx = (b-a)/(n-1);
  hy = (d-c)/(m-1);
  x  = linspace(a,b,n);
  y  = linspace(c,d,m);

  [bx,by,gamma,f,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = local_validacao1(x,y,n,m);

  [ai,bi,ci,di,ei] = coeficientes(hx,hy,kappa,bx,by,gamma,n,m);

  [ai,bi,ci,di,ei,f] = condicoes_contorno(ai,bi,ci,di,ei,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,hx,hy,kappa);
  
  [A,f] = sistema_linear(ai,bi,ci,di,ei,f,n,m);
  
  tDireto = tic();
  u = A\f;
  toc(tDireto);
  
  tSOR = tic();
  [u,iter,er] = sor(ai,bi,ci,di,ei,f,n,m,w,0.00001,1000);
  toc(tSOR);

  grafico_solucao(u,x,y,n,m)

endfunction

function [bx,by,gamma,f,kappa,left,gleft,right,gright,bottom,gbottom,top,gtop] = local_validacao1(x,y,n,m)

  bx = repmat(0.0,1,n*m);
  by = repmat(0.0,1,n*m);
  gamma = repmat(0.0,1,n*m);
  f = repmat(0.0,1,n*m);

  kappa = 1;

  left = 1;
  right = 1;
  bottom = 1;
  top = 1;
  
  gleft = repmat(10,n*m,1);
  gright = repmat(10,n*m,1);
  gbottom = repmat(10,n*m,1);
  gtop = repmat(10,n*m,1);

endfunction