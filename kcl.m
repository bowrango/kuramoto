
rng(0)
clear all
close all

n = 2;

% R = 1e3*ones(n);
R = 1e3*rand(n);
R = (R+R.')/2;
R(1:n+1:end) = 0;

ind = 33e-6;
cap = 750e-12;
f0 = 1/(2*pi*sqrt(ind*cap)); % ~1MHz
omega0 = 2*pi*f0;
T0 = 1/f0;

% opposite polarity coupling
s = (-1).^(0:n-1);
S = diag(s);

Gij = 1 ./ R;
Gij(1:n+1:end) = 0; % handle nans
g = sum(Gij,2);
G = diag(g) - S.*Gij.*S.';

Rcomp = 516.5; % play value
Gcomp = (1/(Rcomp*cap))*eye(n);
% Gcomp = 0;

% 1. KCL Model
A = [zeros(n) eye(n); -(1/(ind*cap))*eye(n) -(1/cap)*G+Gcomp];
B = zeros(size(A,1));
C = eye(2*n);
sys = ss(A, B, C, 0);

Vdd = 3.3; % supply voltage
Vm = Vdd/2;

% TODO Check these units
theta0 = 2*pi*rand(1,n);
v0 = cos(theta0).';
dv0 = -omega0*sin(theta0).';
X0 = [v0; dv0];

P = 20;
dt = T0/P;
t = 0:dt:P*T0;

% U = zeros(length(t), 2*n);
% X = lsim(sys, U, t, X0);

X = initial(sys, X0, t);
v = X(:,1:n);
dv = X(:,n+1:end);

V1 = Vm*v + Vm;

theta1 = unwrap(atan2(-dv/omega0, v), [], 1);

% 2. Kuramoto Model
% K = -G;
% K = (S.*Gij)/(cap*omega0);
K = (S.*Gij.*S.')/(cap*omega0);
K(1:n+1:end) = 0; 
% K = Gij / (cap*omega0);

omegan = sqrt(max(omega0^2 - (g./(2*cap)).^2, 0)); % clip NaNs

kura  = @(t,th) omega0.*ones(n,1) + sum(K.*sin(th'-th), 2);

[~, theta2] = ode45(kura, t, theta0);
V2 = Vm*cos(theta2) + Vm;

figure
subplot(3,1,1)
grid on
plot(t, mod(theta1, 2*pi), 'r'), hold on % KCL
plot(t, mod(theta2, 2*pi), 'b') % Kuramoto
xlabel('time (s)')
ylabel('\theta (rad)')
hold off

subplot(3,1,2)
grid on
plot(t, V1, 'r-o'), hold on
plot(t, V2, 'b-x')
xlabel('time (s)')
ylabel('voltage')

subplot(3,1,3)
grid on
plot(t, abs(V1-V2), 'k'), hold on
xlabel('time (s)')
ylabel('l2')