%es1_6
close all
clear all
clc

a=1; %raggio dell'elica
b=-0.1; %costante dei passi
t=linspace(0,10*pi,500); 
x = a*cos(t); 
y = a*sin(t);
z = b*t;

figure
plot3(x,y,z)

t=linspace(0,20*pi,500);
a=1;
b=0.1;
z=b*t;
x=a*cos(t);
y=a*sin(t);

figure
plot3(x,y,z)