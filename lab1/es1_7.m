%es1_7
close all
clear all
clc

x=linspace(0,1);
y=linspace(0,1);
[X Y]=meshgrid(x,y);
f1= X.*(X-1).*Y.*(Y-1);
f2= X.*(X-1).*sin(X.*8).*Y.*(Y-1).*cos(Y.*8);

figure(11)
mesh(X,Y,f1)
figure(12)
mesh(X,Y,f2)
figure(21)
surf(X,Y,f1)
figure(22)
surf(X,Y,f2)
