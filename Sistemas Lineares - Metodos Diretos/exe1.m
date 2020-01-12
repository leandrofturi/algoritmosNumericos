load bcsstk01.mat
M1 = Problem.A

load 685_bus.mat
M2 = Problem.A

load bcsstk16.mat
M3 = Problem.A

load onetone1.mat
M4 = Problem.A

%load pre2.mat
%M5 = Problem.A

[L1, U1, P1] = lu(M1)
[L2, U2, P2] = lu(M2)
[L3, U3, P3] = lu(M3)
[L4, U4, P4] = lu(M4)
%[L5, U5, P5] = lu(M5)

esparsidade = [nnz(M1), nnz(L1), nnz(P1);
               nnz(M2), nnz(L2), nnz(P2);
               nnz(M3), nnz(L3), nnz(P3);
               nnz(M4), nnz(L4), nnz(P4)]%;
               %nnz(M5), nnz(L5), nnz(P5)]
               

b1 = M1*ones(size(M1)(1,1), 1)
b2 = M2*ones(size(M2)(1,1), 1)
b3 = M3*ones(size(M3)(1,1), 1)
b4 = M4*ones(size(M4)(1,1), 1)
%b5 = M5*ones(size(M5)(1,1), 1)

x1 = M1\b1
x2 = M2\b2
x3 = M3\b3
x4 = M4\b4
%x5 = M5\b5

residuo = [norm(b1-M1*x1, inf);
           norm(b2-M2*x2, inf);
           norm(b3-M3*x3, inf);
           norm(b4-M4*x4, inf)]%;
           %norm(b5-M5*x5, inf)]

erro = [norm(ones(size(M1)(1,1), 1)-x1, inf);
        norm(ones(size(M2)(1,1), 1)-x2, inf);
        norm(ones(size(M3)(1,1), 1)-x3, inf);
        norm(ones(size(M4)(1,1), 1)-x4, inf)]%;
        %norm(ones(size(M5)(1,1), 1)-x5, inf)]