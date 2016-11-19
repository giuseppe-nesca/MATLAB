function [x,iter] = gauss_seidel(A,b,itermax,toll,x0)
% [x,iter] = gauss_seidel(A,b,itermax,toll,x0)
% Input:
%   A la matrice del sistema lineare
%   b vettore termine noto
%   itermax numero massimo iterazioni
%   toll tolleranza richiesta
%   x0 vettore soluzione iniziale
% Output:
%   x soluzione del sistema lineare
%   iter numero di iterazione raggiunte

% Ax=b -> (D+C)x=b -> Dx=b-Cx
D=tril(A); % Gauss-Seidel: triangolare inferiore di A
C=A-D;

% x = invD * (b-Cx)
%   = invD*b - invD*C*x
%   = invD*b - invD*C*(invD*b - invD*C*x)
%   = ... + (-invD*C)^2*x
%   = ... + (-invD*C)^3*x
% B = -invD*C = -invD*(A-D)
%   = -invD*A + I
n=length(b);
B=-inv(D)*A + eye(n); % matrice di iterazione
rho=max(abs(eig(B))); % raggio spettrale
disp(rho)
for iter=1:itermax
    x=D\(b-C*x0);
    if norm(x-x0,'inf')<=toll*norm(x,'inf')
        break % esce dal ciclo for
    end
    x0=x;
end


end