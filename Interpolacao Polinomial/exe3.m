fileID = fopen('exe.txt','w');
fprintf(fileID, 'Q01\n');
for i = 2:7
  x = linspace (-1, 1, i);
  y = zeros(1, i);
  for j = 1:length(x)
    y(j) = 1/(x(j)+10) + x(j)^3 + x(j)^2 - 3;
  endfor
  p = polyfit(x, y, i-1);
  fprintf(fileID, polyout(p, "x"));
  fprintf(fileID, '\n');
  x = linspace (-1, 1, 100);
  if i == 2
    y2 = polyval(p, x);
  endif
  if i == 3
    y3 = polyval(p, x);
  endif
  if i == 4
    y4 = polyval(p, x);
  endif
  if i == 5
    y5 = polyval(p, x);
  endif
  if i == 6
    y6 = polyval(p, x);
  endif
  if i == 7
    y7 = polyval(p, x);
  endif
endfor
xp = linspace (-1, 1, 100);
y = zeros(1, 100);
for j = 1:100
  y(j) = 1/(xp(j)+10) + xp(j)^3 + xp(j)^2 - 3;
endfor
figure();
plot(xp, y, x, y2, x, y3, x, y4, x, y5, x, y6, x, y7);
xlabel ("x");
ylabel ("y");
legend("f(x)", "n = 1", "n = 2", "n = 3", "n = 4", "n = 5", "n = 6");

fprintf(fileID, '\nQ02\n');
x = 6:16;
y = [0.029, 0.052, 0.079, 0.125, 0.181, 0.261, 0.425, 0.738, 1.130, 1.882, 2.812];
figure ();
plot(x, y);
xlabel ("x");
ylabel ("y");
x4 = [x(1), x(4), x(7), x(9), x(11)];
y4 = [y(1), y(4), y(7), y(9), y(11)];
p4 = polyfit(x4, y4, 4);
fprintf(fileID, polyout(p4, "x"));
fprintf(fileID, '\n');
x6 = [x(1), x(3), x(5), x(7), x(9), x(10), x(11)];
y6 = [y(1), y(3), y(5), y(7), y(9), y(10), y(11)];
p6 = polyfit(x6, y6, 6);
fprintf(fileID, polyout(p6, "x"));
fprintf(fileID, '\n');
x8 = [x(1), x(3), x(5), x(6), x(7), x(8), x(9), x(10), x(11)];
y8 = [y(1), y(3), y(5), y(6), y(7), y(8), y(9), y(10), y(11)];
p8 = polyfit(x8, y8, 8);
fprintf(fileID, polyout(p8, "x"));
fprintf(fileID, '\n');
x = linspace (6, 16, 100);
y4 = polyval(p4, x);
y6 = polyval(p6, x);
y8 = polyval(p8, x);
figure ();
plot(x, y4, x, y6, x, y8);
xlabel ("x");
ylabel ("y");
legend("n = 4", "n = 6", "n = 8");
fprintf(fileID, 'f(8.5):\n');
fprintf(fileID, '%f\n', polyval(p4, 8.5));
fprintf(fileID, '%f\n', polyval(p6, 8.5));
fprintf(fileID, '%f\n', polyval(p8, 8.5));
y6 = [x(4), x(5), x(6), x(7), x(8), x(9), x(10)];
x6 = [y(4), y(5), y(6), y(7), y(8), y(9), y(10)];
p6 = polyfit(x6, y6, 6);
fprintf(fileID, 'polinomio inverso: ');
fprintf(fileID, polyout(p6, "x"));
fprintf(fileID, '\n');
fprintf(fileID, 'dia calculado: ');
fprintf(fileID, '%f\n', polyval(p6, 0.3));

fprintf(fileID, '\nQ03\n');
x = linspace (0, 10, 10);
y = sin(2*pi*x/5);
xf = linspace (0, 10, 100);
yfl = interp1(x, y, xf, 'linear');
yfc = interp1(x, y, xf, 'cubic');
yfs = interp1(x, y, xf, 'spline');
xp = linspace (0, 10, 100);
yp = sin(2*pi*xp/5);
figure ();
plot(xp, yp, xf, yfl, xf, yfc, xf, yfs);
xlabel ("x");
ylabel ("sin(2*pi*x/5)");
legend("f(x)", "linear", "cubic", "spline");

fprintf(fileID, '\nQ04\n');
x3 = linspace (-5, 5, 4);
y3 = zeros(1, 4);
for i = 1:4
  y3(i) = 1/(1+0.25*x3(i)^2);
endfor
p3 = polyfit(x3, y3, 3);
x5 = linspace (-5, 5, 6);
y5 = zeros(1, 4);
for i = 1:6
  y5(i) = 1/(1+0.25*x5(i)^2);
endfor
p5 = polyfit(x5, y5, 5);
x10 = linspace (-5, 5, 11);
y10 = zeros(1, 11);
for i = 1:11
  y10(i) = 1/(1+0.25*x10(i)^2);
endfor
p10 = polyfit(x10, y10, 10);
fprintf(fileID, polyout(p3, "x"));
fprintf(fileID, '\n');
fprintf(fileID, polyout(p5, "x"));
fprintf(fileID, '\n');
fprintf(fileID, polyout(p10, "x"));
fprintf(fileID, '\n');
x = linspace (-5, 5, 100);
y3 = polyval(p3, x);
y5 = polyval(p5, x);
y10 = polyval(p10, x);
y = zeros(1, 100);
for i = 1:100
  y(i) = 1/(1+0.25*x(i)^2);
endfor
figure ();
plot(x, y, x, y3, x, y5, x, y10);
xlabel ("x");
ylabel ("f(x)");
legend("f(x)", "n = 3", "n = 5", "n = 10");

x = linspace (-5, 5, 10);
y = zeros(1, 10);
for i = 1:10
  y(i) = 1/(1+0.25*x(i)^2);
endfor
xf = linspace (-5, 5, 100);
yfl = interp1(x, y, xf, 'linear');
yfc = interp1(x, y, xf, 'cubic');
yfs = interp1(x, y, xf, 'spline');
xp = linspace (-5, 5, 100);
for i = 1:100
  yp (i) = 1/(1+0.25*xp(i)^2);
endfor
figure ();
plot(xp, yp, xf, yfl, xf, yfc, xf, yfs);
xlabel ("x");
ylabel ("f(x)");
legend("f(x)", "linear", "cubic", "spline");
fclose(fileID);