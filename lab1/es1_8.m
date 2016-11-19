%es1_8
% I seguenti numeri vengono introdotti in un calcolatore nel quale i numeri vengono rappresentati
%in aritmetica floating-point, con base N = 10 e t = 5 cifre riservate alla mantissa (tecnica di
%arrotondamento (ii)):
%a = 1.483593,
%b = 1.484111.
% Utilizzare il comando chop di MATLAB per determinare il risultato s̄ = ā 	 b̄, ove x̄ denota il
% numero di macchina corrispondente ad x e 	 denota l’operazione di macchina corrispondente
% all’operazione aritmetica della sottrazione. Confrontare s̄ con c = a − b e calcolare l’errore
% relativo corrispondente.


close all
clear all
clc

N=10; %BASE
t=5; %cifre riservate per la MANTISSA
a = 1.483593;
b = 1.484111;
s=a-b;

s_ = chop(chop(a,t)-chop(b,t),t);
err=abs((s-s_)/abs(s))

