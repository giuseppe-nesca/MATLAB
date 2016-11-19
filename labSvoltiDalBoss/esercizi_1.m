disp('***********************************************')
disp('*****************esercizio_1_1*****************')
disp('***********************************************')

x=1:-0.1:0; 

x([1 4 3]) % estrae prima, quarta e terza componente di x
x([1:2:7 10])=zeros(1,5) % sostituisce 0 in posizioni dispari tra 1 e 7 e in posizione 10
x([1 2 5])=[0.5*ones(1,2) -0.3] % sostituisce 0.5 in posizione 1 e 2 e -0.3 in posizione 5
y=x(end:-1:1) % salva in y il vettore x dall'ultima alla prima componente

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_2*****************')
disp('***********************************************')

clear all

A=[1:4;5:8;9:12];

size(A) % restituisce dimensioni della matrice A
B=A.*A % moltiplicazione elemento per elemento
%B=A*A % moltiplicazione riga per colonna: ERRORE di dimensioni
B=A'*A % # colonne A trasposto == # righe A --> B matrice 4x4
B=A*A' % # colonne A == # righe A trasposto --> B matrice 3x3
A(1:2,4)     % seleziona righe 1 e 2, colonna 4 di A
A(:,3)       % seleziona la terza colonna di A (tutte le righe)
A(1:2,:)     % seleziona prima e seconda riga di A (tutte le colonne)
A(:,[2 4])   % seleziona seconda e quarta colonna di A (tutte le righe)
A([2 3 3],:) % seleziona seconda, terza e di nuovo terza riga di A (tutte le colonne)
A(3,2)=A(1,1) % sostituisce A(1,1) in posizione (3,2)
A(1:2,4)=zeros(2,1) % sostituisce 0 nella prima e seconda riga di A, quarta colonna
A(2,:)=A(2,:)-A(2,1)/A(1,1)*A(1,:) % sostituisce alla seconda riga di A, la seconda riga meno la prima moltiplicata per A(2,1)/A(1,1)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_3*****************')
disp('***********************************************')

clear all

D = diag(5*ones(10,1));     % matrice diagonale principale
CS = diag(3*ones(9,1),1);   % matrice codiagonale superiore
CI = diag(-1*ones(9,1),-1); % matrice codiagonale inferiore
B = D+CS+CI;
B([5 8],[6 9]) = 2 % sostituisce 2 in posizione (5,6) (5,9) (8,6) e (8,9)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_4*****************')
disp('***********************************************')

clear all

A = [1:4;2:5;3:6;4:7];
lambda = 10;
P1 = eye(4);
P1(2,2) = lambda;
L1 = P1*A
R1 = A*P1
%
P2 = eye(4);
P2([2 4],:) = P2([4 2],:);
L2 = P2*A
R2 = A*P2
%
P3 = eye(4);
P3(4,2) = lambda;
L3 = P3*A
R3 = A*P3
%
L1a = A;
L1a(2,:) = lambda*L1a(2,:)
R1a = A;
R1a(:,2) = lambda*R1a(:,2)
%
L2a = A;
L2a([2 4],:) = L2a([4 2],:)
R2a = A;
R2a(:,[2 4]) = R2a(:,[4 2])
%
L3a = A;
L3a(4,:) = L3a(4,:)+lambda*L3a(2,:)
R3a = A;
R3a(:,2) = R3a(:,2)+lambda*R3a(:,4)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_5*****************')
disp('***********************************************')

clear all

x = linspace(0,pi/4);
y = tan(x);
figure
plot(x,y)
%
x = linspace(pi/4,pi/2-eps);
y = tan(x);
figure
semilogy(x,y)
%
x = linspace(0.1,100,1000);
y = sqrt((100*(1-0.01*x.^2).^2+0.02*x.^2)./((1-x.^2).^2+0.1*x.^2));
figure
subplot(2,2,1)
plot(x,y)
subplot(2,2,2)
semilogx(x,y)
subplot(2,2,3)
semilogy(x,y)
subplot(2,2,4)
loglog(x,y)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_6*****************')
disp('***********************************************')

close all
clear all

t = linspace(0,10*pi,500);
a = 1;
b = -0.1;
x = a*cos(t);
y = a*sin(t);
z = b*t;
figure
plot3(x,y,z)

t = linspace(0,20*pi,500);
a = 1;
b = 0.1;
x = a*cos(t);
y = a*sin(t);
z = b*t;
figure
plot3(x,y,z)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_7*****************')
disp('***********************************************')

close all
clear all

x = linspace(0,1);
y = linspace(0,1);
[X,Y] = meshgrid(x,y);
f1 = X.*(X-1).*Y.*(Y-1);
figure
mesh(f1)
figure
surf(f1)

f2 = X.*(X-1).*sin(8*X).*Y.*(Y-1).*cos(8*Y);
figure
mesh(X,Y,f2)
figure
surf(X,Y,f2)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc

disp('***********************************************')
disp('*****************esercizio_1_8*****************')
disp('***********************************************')

clear all
close all

a = 1.483593;
b = 1.484111;
t = 5;
s = a-b
s_ = chop(chop(a,t)-chop(b,t),t) % sottrazione di macchina tra numeri di macchina
er = abs(s-s_)/abs(s)

%
disp('Premere un tasto per continuare con l''esercizio successivo...')
pause
%

clc 

disp('***********************************************')
disp('*****************esercizio_1_9*****************')
disp('***********************************************')

clear all
close all

format long e

n = 1:16;
x = 10.^(-n);

f1 = (1-cos(x))./x.^2;
% cos(2*alfa) = cos(alfa)^2-sin(alfa)^2
% 2*alfa = x --> alfa = x/2
% ( 1 - cos(x/2).^2 + sin(x/2).^2 )./x.^2
% ( 1 - 1 + sin(x/2).^2 + sin(x/2).^2 )./x.^2
% 2*sin(x/2).^2 ./ x.^2
% 2*( sin(x/2)./ x ).^2
% 1/2*( sin(x/2) ./ (x/2) ).^2
f1_ex = 1/2*(sin(x/2)./(x/2)).^2;
disp([f1' f1_ex'])
er1 = abs(f1-f1_ex)./abs(f1_ex);
disp([x' er1'])
figure
loglog(x,er1)

%
disp('Premere un tasto per continuare con la funzione successiva...')
pause
%

f2 = (1-exp(x))./x;
% svliuppo in serie di Taylor/McLaurin
% exp(x) = 1 + x + x.^2/2 
%            + x.^3/factorial(3) + ...
% ( 1 - 1 - x - x.^2/2 - x.^3/factorial(3) - ...)./x
% (- x - x.^2/2 - x.^3/factorial(3)-..)/x
% - 1 - x/2 - x.^2/facorial(3) - ...
f2_ex = 0;
for i = 1:16
    f2_ex = f2_ex - x.^(i-1)/factorial(i);
end    
disp([f2' f2_ex'])
er2 = abs(f2-f2_ex)./abs(f2_ex);
disp([x' er2'])
figure
loglog(x,er2)

%
disp('Premere un tasto per continuare con la funzione successiva...')
pause
%

f3 = (1-sqrt(1-x.^2));
% razionalizzando
% (1-sqrt(1-x.^2)) .* (1+sqrt(1-x.^2)) ./ (1+sqrt(1-x.^2))
% (1-(1-x.^2)) ./ (1+sqrt(1-x.^2))
% x.^2 ./ (1+sqrt(1-x.^2))
f3_ex = x.^2./(1+sqrt(1-x.^2));
disp([f3' f3_ex'])
er3 = abs(f3-f3_ex)./abs(f3_ex);
disp([x' er3'])
figure
loglog(x,er3)

%
disp('Premere un tasto per continuare con la funzione successiva...')
pause
%

f4 = ((x+1).^2-1)./x;
% svolgendo il quadrato
% (x.^2+2*x+1 - 1) ./ x
% (x.^2+2*x) ./ x
% x+2
f4_ex = x+2;
disp([f4' f4_ex'])
er4 = abs(f4-f4_ex)./abs(f4_ex);
disp([x' er4'])
figure
loglog(x,er4)


disp('Fine del foglio Esercitazione 1')