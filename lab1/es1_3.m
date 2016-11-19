%es1_3
clear all
clc

D=diag(5*ones(10,1));
DS=diag(3*ones(9,1),1);
DI=diag(-1*ones(9,1),-1);
B=D+DS+DI;
B([5 8],[6 9])=2;
disp(B)