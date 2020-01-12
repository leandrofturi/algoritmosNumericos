function [error_x, Keps, Kerror_b, nresiduo] = avaliacond(n);
  H = hilb(n);
  b = H*ones(n, 1);
  x = H\b;
  dx = ones(n, 1)-x;
  db = H*dx;
  nb = norm(db,inf)/norm(b,inf);
  error_x = norm(dx,inf)/norm(x,inf);
  K = cond(H);
  Keps = K*eps;
  Kerror_b = K*nb;
  nresiduo = norm(b-H*x,inf);
endfunction