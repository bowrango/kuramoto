
% "Inverse optimal annealing" viewpoint for SDRE-controlled gradient flow

clear
close all

rng default

n = 8;

J = randn(n);
J = (J+J.')/2;
J(1:n+1:end) = 0;

qb = qubo(J);
gt = solve(qb);

Qe = eye(n);
R  = 1e8*eye(n);

theta0 = 2*pi*rand(n,1) - pi;

params.J  = J;
params.Qe = Qe;
params.R  = R;

tspan = [0 10];
rhs = @(t,x) ising_sdre_rhs(t, x, params);
[t, X] = ode45(rhs, tspan, theta0);

N = length(t);
H = zeros(N,1);
Enorm = zeros(N,1);
u_hist = zeros(N, n);
B = eye(n);
for k = 1:N
    xk = X(k,:).';

    [H(k), ek, Mk] = ising_energy_and_gradient(xk, J);

    A = -Mk;
    Qx = Mk.'*Qe*Mk + 1e-6*eye(n);
    P = care(A, B, Qx, R);
    Kx = R \ (B.' * P);

    u_hist(k,:) = (-Kx * xk).';
    Enorm(k) = norm(ek, 2);
end

dt = [diff(t); t(end) - t(end-1)];
instEn = sum(u_hist.^2, 2);
cumEn = cumsum(instEn .* dt);

figure;
subplot(4,1,1);
plot(t, X, 'LineWidth', 1.0);
ylabel('\theta_i(t)');
title('Phases under SDRE gradient flow');
grid on;

subplot(4,1,2);
hold on
plot(t, H, 'LineWidth', 1.5);
yline(gt.BestFunctionValue, '--', LineWidth=1.5);
hold off
ylabel('H(\theta(t))');
title('Energy H(\theta)');
grid on;

subplot(4,1,3);
plot(t, Enorm, 'LineWidth', 1.5);
ylabel('||e(\theta)||_2');
title('Residual ||e(\theta)||');
grid on;

subplot(4,1,4);
plot(t, instEn, 'LineWidth', 1.5);
xlabel('time t');
ylabel('||u(t)||^2');
title('Control Effort');
grid on;

function dtheta = ising_sdre_rhs(~, theta, params)

J  = params.J;
Qe = params.Qe;
R  = params.R;

n = length(theta);

[~, e, M] = ising_energy_and_gradient(theta, J);

A = -M;
B = eye(n);

% Map residual cost into an effective quadratic cost on theta
Qx = M.'*Qe*M + 1e-6*eye(n);

P  = care(A, B, Qx, R);
Kx = R \ (B.' * P);

u = -Kx * theta;
dtheta = -e + u;
end

function [H, e, M] = ising_energy_and_gradient(theta, J)

z = cos(theta);

H = -0.5 * (z.' * J * z);

c  = z;
s  = sin(theta);
Jc = J * c;
e  = s .* Jc;

n = length(theta);
M = -(s * s.') .* J;

diagVals = c .* Jc;
M(1:n+1:end) = diagVals;
end