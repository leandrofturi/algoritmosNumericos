fileID = fopen('exe.txt','w');
figure ();
hold on;

V = [2.6500   2.6500   2.7000   2.7000   2.7500   2.7500   2.8500   2.8500   2.9000   2.9000   2.9500   2.9500   3.0000   3.0000]
T = [6.8500   6.8000   6.7000   6.3000   6.3300   6.2000   5.9000   5.8200   5.8000   5.8000   6.1500   6.0000   6.3000   6.1500]

plot(V, T, '*');
colors = ['r', 'g', 'y', 'm', 'c', 'k'];
for(i = 1:6)
  fprintf(fileID, '\n');
  p = polyfit(V, T, i);
  fprintf(fileID, polyout(p, "x"));
  r = (T - polyval(p, V));
  Sr = dot(r, r);
  r2 = 1 - (Sr/(dot(T, T) - ((sum(T)^2)/14)));
  phi2 = Sr / (14-i-1);
  fprintf(fileID, '\n');
  fprintf(fileID, "r2 = %f", r2);
  fprintf(fileID, '\n');
  fprintf(fileID, "phi2 = %f", phi2);
  fprintf(fileID, '\n');
  fprintf(fileID, "T(2.8) = %f", polyval(p, 2.8));
  
  x = linspace(V(1), V(14), 100);
  y = polyval(p, x);
  plot(x, y, colors(i));
  fprintf(fileID, '\n');
endfor

xlabel ("V");
ylabel ("T");
legend("T(V)", "1", "2", "3", "4", "5", "6");
hold off;

S = [1.3  1.8  3  4.5  6  8  9];
V = [0.07  0.13  0.22  0.275  0.335  0.35  0.36];

x = 1./S;
y = 1./V;
fprintf(fileID, '\n');
fprintf(fileID, 'Caso 1:\n');
p = polyfit(x, y, 1);
Vm = 1/p(2);
Ks = p(1)*Vm;
fprintf(fileID, 'Vm = %f, Ks = %f', Vm, Ks);
fprintf(fileID, '\n');
fprintf(fileID, 'S(7) = %f', (Vm*7)/(Ks+7));
r = (y - polyval(p, x));
Sr = dot(r, r);
r2 = 1 - (Sr/(dot(y, y) - ((sum(y)^2)/7)));
phi2 = Sr / (7-1-1);
fprintf(fileID, '\n');
fprintf(fileID, "r2 = %f", r2);
fprintf(fileID, '\n');
fprintf(fileID, "phi2 = %f", phi2);
fprintf(fileID, '\n');
x1 = linspace(1.3, 9, 100);
y1 = (Vm*x1);
y1 = y1./(Ks.+x1);
w1 = linspace(1/1.3, 1/9, 100);
z1 = polyval(p, w1);


x = 1./(S.^2);
y = 1./V;
fprintf(fileID, '\n');
fprintf(fileID, 'Caso 2:\n');
p = polyfit(x, y, 1);
Vm = 1/p(2);
Ks = sqrt(p(1)*Vm);
fprintf(fileID, 'Vm = %f, Ks = %f', Vm, Ks);
fprintf(fileID, '\n');
fprintf(fileID, 'S(7) = %f', (Vm*(7^2))/((Ks^2)+(7^2)));
r = (y - polyval(p, x));
Sr = dot(r, r);
r2 = 1 - (Sr/(dot(y, y) - ((sum(y)^2)/7)));
phi2 = Sr / (7-1-1);
fprintf(fileID, '\n');
fprintf(fileID, "r2 = %f", r2);
fprintf(fileID, '\n');
fprintf(fileID, "phi2 = %f", phi2);
fprintf(fileID, '\n');
x2 = linspace(1.3, 9, 100);
y2 = (Vm*(x2.^2));
y2 = y2./((Ks^2).+(x2.^2));
w2 = linspace(1/(1.3^2), 1/(9^2), 100);
z2 = polyval(p, w2);

figure();
plot(S, V, '*', x1, y1, x2, y2)
xlabel ("S");
ylabel ("V");
legend("V(S)", "caso 1", "caso 2");

figure();
plot(1./S, 1./V, '*', w1, z1)
xlabel ("S");
ylabel ("V");
legend("1/S x 1/V", "estimado");

figure();
plot(1./(S.^2), 1./V, '*', w2, z2)
xlabel ("S");
ylabel ("V");
legend("1/S² x 1/V", "estimado");

fclose(fileID);