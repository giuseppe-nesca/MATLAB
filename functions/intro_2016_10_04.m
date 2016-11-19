% --> inserisce commenti
%% variabili scalari
a=3
b=4; % ; non stampa operazione

c=a+b;
disp(c) % stampa valore di c

% case sensitive:
Aa=67;
aA=55;

% precedenze operazioni
d=a+b*c
dd=(a+b)*c
ddd=((a+b)*c)^2

pi %pi greco
pi=1;
clear pi %cancella variabile

format long % 15 cifre decimali
pi
format long e % potenza di 10
disp(0.0000123456)
format short e
disp(0.0000123456)
format short % 4 cifre decimali (default)

l=log(5); % logaritmo naturale (in base e)
l1=log10(1000); % logaritmo in base 10
l2=log2(64); % logaritmo in base 2

s=sqrt(2); % radice quadrata

% altre funzioni predefinite: exp, trigonometriche (sin, cos, tan, ...)

who % lista variabili assegnate
whos a
whos a* % tutte quelle che iniziano per a

clear all % cancella tutte le variabili

clc % cancella linea di comando

%% variabili vettoriali
v=[1 2,3]; %spazio/virgola -> vettore riga
w=[3;4;5]; %; -> vettore colonna
length(v) % lunghezza vettore
vT=v'; % vettore trasposto
w+vT  % dimensioni ok: 3x1 + 3x1
% w*vT  % errore: 3x1 * 3x1
v*w % ok: 1x3 * 3x1 -> 1x1
w.*vT % ok (elemento per elemento)

a=1;
b=10;
z=a:b  % vettore di interi da 1 a 10
zz=a:0.5:b % vettore da 1 a 10 con passo 0.5
zzz=b:-1:a % vettore con passo negativo (da 10 a 1)

y=linspace(a,b,5) % vettore di 5 elementi equispaziati tra a e b
y(1) % primo elemento
y(end) % y(length(y)) % y(N)
y([1 3 5])
y(2:2:end) % componenti in pos pari
y(3)=0 % sostituisce 0 in pos 3
y=y(end:-1:1) % gira y al contrario

u=ones(3,1);  % vettore colonna 3x1 di 1
z=zeros(1,4) % vettore riga 1x4 di 0
r=rand(1,3) % vettore di numeri pseudocasuali


%% matrici
%% matrici
A=[1 2 3; 4 4 5; 9 8 7];
B=rand(3,2)
size(A) % dimensioni matrice

B*A % errore: 3x2 * 3x3
C=A*B % 3x3 * 3x2 --> 3x2
D=C.*B % moltiplicazione elem X elem

A^2 % == A*A (deve essere quadrata)
A.^2 % == A.*A

det(A) % determinante
rank(A) % rango
eig(A) % autovalori
[autovett,autoval]=eig(A)

diag(autoval) % estrae diagonale
diag(v) % mette v sulla diagonale

A(1,1) % elemento riga 1, colonna 1
A(end,end)
r=2;
c=3;
A(r,c) % riga r, colonna c
A(r,:) % tutta la riga r
A(:,c) % tutta la colonna c
A(r,c)=0 % 0 nella riga r, colonna c
% scambia colonne 1 e 2:
A(:,[1 2])=A(:,[2 1])
% eliminare riga r:
A(r,:)=[]
% aggiungere una colonna a B
B=[B ones(3,1)]

%% nel salvataggio NON usare
% + - / * . spazio numeri (da soli)
% usate _ lettere numeri (non soli)

