%% Esercizio_2_1
clear all
clc
disp('***********************************************')
disp('*****************esercizio_2_1*****************')
disp('***********************************************')

format long e
for n = 5:5:15 % [5 10 15] % ciclo su dimensioni matrice
    A = hilb(n); % matrice di Hilbert
    b = sum(A,2); %A*ones(n,1); % somma elementi su righe
    x = A\b; % \ risolve il sistema lineare
    x_ex = ones(n,1);
    % err relativo in norma euclidea (2):
    err_rel = norm(x-x_ex,2)/norm(x_ex,2);
    disp([n err_rel cond(A)]) % cond dà numero di condizionamento
end

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

%% Esercizio_2_2
clc
format short
disp('***********************************************')
disp('*****************esercizio_2_2*****************')
disp('***********************************************')

clear all

A = [11 2 3 4; 0 8 2 3; 0 0 4 0;1 0 0 5];
b = sum(A,2); % solizione vettore unitario

[L,U,P] = lu(A); %P*A=L*U

% A x = b
% P A x = P b -> L U x = P b
% L y = P b  ,  U x = y
y = L\(P*b);
x = U\y

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

%% Esercizio_2_3
clc
disp('***********************************************')
disp('*****************esercizio_2_3*****************')
disp('***********************************************')

clear all

A=rand(3);
if det(A)==0 % matrice singolare
    disp('ERRORE: A singolare')
    return % termina il programma
end

[L,U,P]=lu(A); % P*A=L*U

% A = invP * L * U
% invA = invU * invL * P
invA = inv(U)*inv(L)*P
inv(A)

max(max(invA-inv(A))) % un 'max' per righe, uno per colonne

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

%% Esercizio_2_4
clc
disp('***********************************************')
disp('*****************esercizio_2_4*****************')
disp('***********************************************')

clear all

B = rand(5);
if det(B)==0
    disp('ERRORE: B singolare')
    return
end
A=B'*B;
if isequal(A,A')
    disp('A simmetrica')
else
    disp('ERRORE: A non simmetrica')
    return
end
eig(A) % autoval > 0 allora A è def pos

R = chol(A); % decomposizione Choleski

% A = R'*R
% invA = inv(R'*R)
%      = inv(R) * inv(R')
%      = inv(R) * inv(R)'
invR=inv(R);
invA=invR*invR'
inv(A)

b=sum(A,2); % soluz vettore unitario
% A x = b
% R'*R*x = b
% R'*y=b  ,  R*x=y
y=R'\b;
x=R\y


%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

%% Esercizio_2_5
clc
disp('***********************************************')
disp('*****************esercizio_2_5*****************')
disp('***********************************************')

clear all

A = rand(6);
b = sum(A,2);
[L,U,P] = lu(A); %P*A=L*U

% iniziazzazione matrice di soluzioni in ogni colonna salviamo x_i
x = zeros(6,4);
% Ax=b -> PAx=Pb -> LUx=Pb
% Ly=Pb  ,  Ux=y
y=L\(P*b);
x(:,1)=U\y;
for i = 2:4
    b = b/i;
    y = L\(P*b);
    x(:,i) = U\y;
end 
disp(x)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

%% Esercizio_2_6
clc
disp('***********************************************')
disp('*****************esercizio_2_6*****************')
disp('***********************************************')

clear all

A = rand(5);
N=4; % esponente: A^N * z = b
b = sum(A^N,2);

[L,U,P] = lu(A); % P*A=L*U

% A^4 z = b
% A*A*A*A*z = b
% A*z=y , A*A*A*y=b
% A*z=y , A*y=w , A*A*w=b
% A*z=y , A*y=w , A*w=x , A*x=b

tn = P*b;
for i = 1:N
    v = L\(P*b);
    b = U\v; % nuovo termine noto
end
z=b; % l'ultimo termine noto è la soluzione
disp(z)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

%% Esercizio_2_7
clc
disp('***********************************************')
disp('*****************esercizio_2_7*****************')
disp('***********************************************')

clear all

itermax = 100;
toll = 1e-7; % 1*10^-7

A1 = [1 -2 2; -1 1 -1; -2 -2 1];
b1 = sum(A1,2);
n = length(b1);
x0 = zeros(n,1);
disp('Jacobi:')
tic % fa partire un contatore del tempo
[xJ1,iterJ1] = jacobi(A1,b1,itermax,toll,x0)
toc % termina il contatore del tempo
disp('Gauss-Seidel:')
tic % fa partire un contatore del tempo
[xGS1,iterGS1] = gauss_seidel(A1,b1,itermax,toll,x0)
toc % termina il contatore del tempo

%
disp('Premere un tasto per continuare con il sistema successivo...')
pause
%

A2 = [4 0 2/5; 0 5 2/5; 5/2 2 1];
b2 = sum(A2,2);
n = length(b2);
x0 = zeros(n,1);
disp('Jacobi:')
tic % fa partire un contatore del tempo
[xJ2,iterJ2] = jacobi(A2,b2,itermax,toll,x0)
toc % termina il contatore del tempo
disp('Gauss-Seidel:')
tic % fa partire un contatore del tempo
[xGS2,iterGS2] = gauss_seidel(A2,b2,itermax,toll,x0)
toc % termina il contatore del tempo

%
disp('Premere un tasto per continuare con il sistema successivo...')
pause
%

A3 = [2 -1 1; 2 2 2; -1 -1 2];
b3 = sum(A3,2);
n = length(b3);
x0 = zeros(n,1);
disp('Jacobi:')
tic % fa partire un contatore del tempo
[xJ3,iterJ3] = jacobi(A3,b3,itermax,toll,x0)
toc % termina il contatore del tempo
disp('Gauss-Seidel:')
tic % fa partire un contatore del tempo
[xGS3,iterGS3] = gauss_seidel(A3,b3,itermax,toll,x0)
toc % termina il contatore del tempo

%
disp('Premere un tasto per continuare con il sistema successivo...')
pause
%

A = [3 -1 0;-1 3 -1;0 -1 3];
B = [0 0 -1;0 -1 0;-1 0 0];
A4 = [A B;B A];
b4 = sum(A4,2);
n = length(b4);
x0 = zeros(n,1);
disp('Jacobi:')
tic % fa partire un contatore del tempo
[xJ4,iterJ4] = jacobi(A4,b4,itermax,toll,x0)
toc % termina il contatore del tempo
disp('Gauss-Seidel:')
tic % fa partire un contatore del tempo
[xGS4,iterGS4] = gauss_seidel(A4,b4,itermax,toll,x0)
toc % termina il contatore del tempo

    