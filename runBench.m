
rng default

K = 1;
Ks = 0.1;
Kn = 0.01;

N = 2.^(1:5);
numTrials = 100;
colors = rand([1,3,length(N)]);

for n = N
for jj = 1:numTrials
A = rand(n);
A(1:n+1:end) = 0;
A = A+A.';
g = graph(A);
J = -A;

[cuts, steps] = run(J, K, Ks, Kn);
semilogy(steps, cuts, Color=colors(:,:,find(n==N)))
hold on
grid on
end
n
end

function [cuts, T] = run(J, K, Ks, Kn)
% Run MaxCut simulation

nOsc = size(J,1);

x0 = rand(nOsc,1);

tstop = 5;
dt = 1e-1;

a1.k = (K-1)/tstop;
coup = @(t, args) 1 + t*args.k;

drift = @(t,X) phaseModel(X, coup(t, a1), Ks, J);
diffusion = @(t,X) Kn*eye(nOsc);

model = sde(drift, diffusion, StartState=x0);
[X, T] = simulate(model, round(tstop/dt), DeltaTime=dt);

cuts = zeros(size(T));
for k = 1:length(T)
    x = X(k,:);
    % Ising function
    % Convert phases to spins
    s = ones(nOsc,1);
    odds = find(mod(round(x), 2));
    s(odds) = -1;
    % Convert to binary
    bits = (s+1)/2;
    % Compute MaxCut. Using spins its (s.'*J*s)/2
    cuts(k) = (~bits.'*(-J)*bits);
end
end

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

