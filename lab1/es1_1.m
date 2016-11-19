%es1_1
clear all
%clc

x=([1:-0.1:0]) %vettore riga da 1 a 0 procedendo per x-0.1 alla volta

x([1 4 3]); %colonne 1 4 3 di x

%x=([1:2:7 10])=zeros(1,5)

x([1 2 5])=[0.5*ones(1,2) -0.3]; %colonne 1 2 5 di x uguali a 0.5* vettore di uni dim1,2 e -0.3

y=x(end:-1:1) %y= x letto da destra a sinistra