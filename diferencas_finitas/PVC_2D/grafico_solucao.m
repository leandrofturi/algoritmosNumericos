function grafico_solucao(u,x,y,n,m)

  figure();
  z = zeros(m, n);
  I = 1;
    for j = 1:m
      for i = 1:n
        z(j, i) = u(I);
        I = I+1;
      end;
    end;
  mesh(x,y,z);
  xlabel('x');
  ylabel('y');
  zlabel('z');
  
endfunction