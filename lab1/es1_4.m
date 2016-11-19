%es1_4
clear all
clc

A=[1:4;2:5;3:6;4:7];
lambda=10;
P1=diag(ones(4,1));
P1(2,2)=lambda;
P2=eye(4);
P3=eye(4);
P3(4,2)=lambda;

L1=P1*A;    %A(2,:) è moltiplicato per 10
R1=A*P1;    %A(:,2) è moltiplicato per 10
L11=A;
L11(2,:)=A(2,:).*lambda;
R11=A;
R11(:,2)=A(:,2).*lambda;

L2=P2*A;    %L2==A
R2=A*P2;    %R2==A

L3=P3*A;
R3=A*P3;
L3=A;
L3(4,:)=A(4,:)+(A(2,:)*lambda);
R3=A;
R3(:,2)=A(:,2)+(A(:,4)*lambda);




