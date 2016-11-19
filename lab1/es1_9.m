%es1_9

%  Valutare le espressioni
% 1 − cos(x)
% f1 (x) =
%  ,
% x2
% 1 − e
% x
% f2 (x) =
%  ,
% xp
% f3 (x) = 1 − 1 − x2 ,
% (x + 1)2 − 1
% f4 (x) =
% x
% in x = 10−n per n = 1, 2, ..., 16. Successivamente riformulare le funzioni assegnate al fine di
% evitare il fenomeno della cancellazione numerica e, assumendo come valori esatti quelli che
% si ottengono mediante la riformulazione proposta, calcolare i corrispondenti errori relativi e
% confrontarli con la precisione di macchina. Stampare e rappresentare graficamente per ogni
% valore di x l’errore relativo corrispondente.

close all
clear all
clc

n=linspace(1,16);
x=10.^n;
f1=(1-cos(x))./x.^2;
f2=(1-exp(x))./x;
f3=1-sqrt(1-x.^2);
f4=((x+1).^2-1)./x;

