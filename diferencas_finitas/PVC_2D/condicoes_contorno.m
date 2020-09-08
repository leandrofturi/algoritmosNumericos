function [a,b,c,d,e,f] = condicoes_contorno(a,b,c,d,e,f,n,m,left,gleft,right,gright,bottom,gbottom,top,gtop,hx,hy,kappa)

  if(((left != 1) && (left != 2) && (left != 3)) ||
     ((right != 1) && (right != 2) && (right != 3)) ||
     ((bottom != 1) && (bottom != 2) && (bottom != 3)) ||
     ((top != 1) && (top != 2) && (top != 3)))
    printf("Erro na Condicao de contorno");
  endif

  % Condicoes de contorno para fluxo prescrito e condicao mista
  % Condicoes de contorno left
  if(left == 2)
    for I = (0:m-1)*n+1
      a(I) = a(I) + b(I);
      f(I) = f(I) + b(I)*hx*gleft(I,1)/kappa;
      b(I) = 0;
    endfor
  elseif(left == 3)
    for I = (0:m-1)*n+1
      a(I) = a(I) + b(I)*(1-(hx*gright(I,2)/gright(I,1)));
      f(I) = f(I) + b(I)*hx*gleft(I,3)/gleft(I,1);
      b(I) = 0;
    endfor
  endif
  
  % Condicoes de contorno right
  if(right == 2)
    for I = (1:m)*n
      a(I) = a(I) + c(I);
      f(I) = f(I) + c(I)*hx*gright(I,1)/kappa;
      c(I) = 0;
    endfor
  elseif(right == 3)
    for I = (1:m)*n
      a(I) = a(I) + c(I)*(1-(hx*gright(I,2)/gright(I,1)));
      f(I) = f(I) - c(I)*hx*gright(I,3)/gright(I,1);
      c(I) = 0;
    endfor
  endif

  % Condicoes de contorno bottom
  if(bottom == 2)
    for I = 1:n
      a(I) = a(I) + d(I);
      f(I) = f(I) + d(I)*hy*gbottom(I,1)/kappa;
      d(I) = 0;
    endfor
  elseif(bottom == 3)
    for I = 1:n
      a(I) = a(I) + d(I)*(1-(hy*gbottom(I,2)/gbottom(I,1)));
      f(I) = f(I) - d(I)*hy*gbottom(I,3)/gbottom(I,1);
      d(I) = 0;
    endfor
  endif

  % Condicoes de contorno top
  if(top == 2)
    for I = ((m-1)*n+1):(m*n)
      a(I) = a(I) + e(I);
      f(I) = f(I) + e(I)*hy*gtop(I,1)/kappa;
      e(I) = 0;
    endfor
  elseif(top == 3)
    for I = ((m-1)*n+1):(m*n)
      a(I) = a(I) + e(I)*(1-(hy*gtop(I,2)/gtop(I,1)));
      f(I) = f(I) + e(I)*hy*gtop(I,3)/gtop(I,1);
      e(I) = 0;
    endfor
  endif

  % Condicoes de contorno para valor prescrito
  % Condicoes de contorno left
  if(left == 1)
    for I = (0:m-1)*n+1
      b(I) = c(I) = d(I) = e(I) = 0;
      a(I) = 1;
      f(I) = gleft(I,1);
    endfor
  endif

  % Condicoes de contorno right
  if(right == 1)
    for I = (1:m)*n
      b(I) = c(I) = d(I) = e(I) = 0;
      a(I) = 1;
      f(I) = gright(I,1);
    endfor
  endif

  % Condicoes de contorno bottom
  if(bottom == 1)
    for I = 1:n
      b(I) = c(I) = d(I) = e(I) = 0;
      a(I) = 1;
      f(I) = gbottom(I,1);
    endfor
  endif

  % Condicoes de contorno top
  if(top == 1)
    for I = ((m-1)*n+1):(m*n)
      b(I) = c(I) = d(I) = e(I) = 0;
      a(I) = 1;
      f(I) = gtop(I,1);
    endfor
  endif

endfunction