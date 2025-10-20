function dxdt = phaseModel(x, K, Ks, J)
% Kuramoto (4.6)
% The phases x are assumed to be relative to a reference oscillator
n = length(x);
dxdt = zeros(n,1);
for ii = 1:n
    dxdt(ii) = -K*( J(ii,:)*sin(pi*(x(ii) - x)) );
end
dxdt = (dxdt - Ks*sin(2*pi*x))/pi;
end