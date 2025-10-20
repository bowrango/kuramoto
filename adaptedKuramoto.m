
clear all
close all

nOsc = 3;

% Problem matrix
rng default
Q = rand(nOsc);
Q = Q+Q.';
Q(1:nOsc+1:end) = 0;

% Classical tabu search
g = graph(Q);
qb = maxcut2qubo(g);
gt = solve(qb);
gt.BestX
gt.BestFunctionValue

% Ising matrix
J = -Q;

tstop = 10;
dt = 1e-2;

% Coupling schedule (ramp)
% The true K of each oscillator depends on PPV and perturbation amplitude
% from other oscillators
K = 1;
a1.k = (K-1)/tstop;
coupling = @(t, args) 1 + t*args.k;

% Sync schedule (square wave)
a2.T = tstop/20;
sync = @(t, args) 1+2*tanh(10*cos(2*pi*t/args.T));

drift = @(t,X) phaseModel(X, coupling(t, a1), sync(t, a2), J);

% Noise schedule (constant)
% Curious how this relates to Boltzmann re: Eq 4.17
Kn = 0.1;
diffusion = @(t,X) Kn*eye(nOsc);


x0 = pi*rand(nOsc,1);
model = sde(drift, diffusion, StartState=x0);

[X, T] = simulate(model, tstop/dt, DeltaTime=dt);

% Compute MaxCut value at each time step
% Phases 0/pi -> +1/-1
cuts = zeros(size(T));
for k = 1:length(T)
    mask = true(nOsc,1);
    odds = find(mod(round(X(k,:)), 2));
    mask(odds) = false;
    cuts(k) = -(~mask.'*J*mask);
end

tiledlayout(3,1)

nexttile
plot(T, X)
grid on
ylabel('phases (\pi)')

nexttile
hold on
grid on
yline(-gt.BestFunctionValue, LineWidth=2)
plot(T, cuts)
ylabel('cut value')
hold off

nexttile
hold on
grid on
plot(T, coupling(T, a1))
plot(T, sync(T, a2))
xlabel('time (cycles)');
hold off

function dxdt = phaseModel(x, K, Ks, J)
% Adapted Kuramoto

% FIXME Lyapunov analysis
% tanh(sin()) used for coupling changes the cos() term in (4.7) to
% triangle function (see page 77). Should be dEdt <= 0
shift = x - x.';
E = -K*sum(J(:).*cos(shift(:))) + Ks*sum(cos(2*x));

n = length(x);
dxdt = zeros(n,1);
for ii = 1:n
    % Coupling, sync, and normalize
    % 4.16
    dxdt(ii) = (-K*J(ii,:)*tanh(10*sin(pi*(x(ii)-x))) - Ks*sin(2*pi*x(ii)))/pi;

    % Basic Kuramoto (4.3 and 4.4)
    % dxdt(ii) = -K*J(ii, :)*sin(x(ii) - x);
    % E(ii) = -K*J(ii, :)*cos(x(ii) - x);
end
end
