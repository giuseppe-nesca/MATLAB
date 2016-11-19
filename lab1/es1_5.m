%es1_5
clear all
close all
clc
x=linspace(0,pi/4);
y=tan(x);
figure(1)
plot(x,y,'r')

x=linspace(pi/4,pi/2);
y=tan(x);
figure(2)
subplot(2,1,1)
plot(x,y,'g')
subplot(2,1,2)
semilogy(x,y,'k')

x=linspace(0.1,100);
y=((100*(1-0.01*x.^2).^2+0.02*x.^2)./((1-x.^2).^2+0.1*x.^2)).^(1/2);
figure(3)
subplot(2,2,1)
plot(x,y,'r')
subplot(2,2,2)
semilogx(x,y)
subplot(2,2,3)
semilogy(x,y)
subplot(2,2,4)
loglog(x,y)