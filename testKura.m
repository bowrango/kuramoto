nOsc = 2;

% J = [ 0 0.9294 0.1682 0.2574 0.1395 0.3284
% 0.9294 0 0.8639 0.8428 1.0879 0.0015
% 0.1682 0.8639 0 1.2203 1.3125 0.0177
% 0.2574 0.8428 1.2203 0 0.8108 0.5174
% 0.1395 1.0879 1.3125 0.8108 0 1.6862
% 0.3284 0.0015 0.0177 0.5174 1.6862 0];

J = rand(nOsc);
J = (J+J.')/2;
J(1:nOsc+1:end)=0;


% g = graph(A);
% 
% qb = maxcut2qubo(g);
% gt = solve(qb);

% J = -A;
% h = -qb.LinearTerm;

tstop = 5; tstep = 1e-3;
As = 2; Ac = 5; An = 0.1;
% F = @(t,X) KuramotoSin(X, Ac*t/tstop, As, nOsc, J);
F = @(t,X) phaseModel(X, Ac*t/tstop, As, nOsc, J);
G = @(t,X) An*eye(nOsc);

rng(0)
x0 = rand(nOsc, 1);

obj = sde(F, G, 'StartState', x0);
[X, T] = simulate(obj, tstop/tstep, 'DeltaTime', tstep);

E = zeros(size(T));
cuts = zeros(size(T));
ising = zeros(size(T));
for k = 1:length(T)
    x = X(k,:)';
    % Ising function
    % Convert phases to spins
    s = ones(nOsc,1);
    odds = find(mod(round(x), 2));
    s(odds) = -1;
    % Convert to binary
    bits = (s+1)/2;
    % Compute MaxCut. Using spins its (s.'*J*s)/2
    ix = find(mod(round(X(k,:)), 2));
    cuts(k) = -(s.'*J*s)/2;
end

figure; plot(T, X);
xlabel('time'); ylabel('phase (\pi)');
box on; grid on;

% figure
% plot(T, cuts)
% yline(-gt.BestFunctionValue)


function fout = KuramotoSin(x, Ac, As, n, J)
for c = 1:n
fout(c, 1) = Ac*J(c, :) * sin(pi*(x(c) - x));
end
fout = (fout - As * sin(2*pi*x))/pi;
end