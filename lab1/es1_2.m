%es1_2
clear all
clc

A=[1:1:4; 5:1:8; 9:1:12;];
disp(A)
disp(A') %A' == A trasposta

size(A); %printa le dimensioni di A
%B=A*A; %moltiplicazione riga per colonna: ERRORE per dimensioni non idonee
B=A.*A; %prodotto di matrici
B=A'*A; %prodotto A trasposta per A
B=A*A'; %prodotta A per A trasposta
A(1:2,4)%mostra righe 1 e 2 , colonna 4
A(:,3) %mostra tutte righe, colonna 3
A(1:2,:) %mostra righe 1 e 2, tutte le colonne
A(:,[2 4]) %mostra tutte le righe , colonna 2 e 4
A([2 3 3],:) %mostra righe 2, 3, 3 , colonne tutte
