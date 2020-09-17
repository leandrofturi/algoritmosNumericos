function er = pos_processamento(x, s, u);
  figure();
  hold on;
  plot(x, s, 'r');
  plot(x, transpose(u), 'c');
  legend("real", "aproximada");
  hold off;
  er = norm(abs(u-s), inf);
endfunction
