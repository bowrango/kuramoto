
close all

% There can be muliple stable non-global stationary distributions. How do 
% we guide the system with u,K,Ks? Kuramoto as SDE Decompose station distribution

rng(0)

% FIXME control gain
% TODO natural frequencies
% TODO B matrix

% [1] https://arxiv.org/pdf/2504.18399
% [2] https://arxiv.org/pdf/2503.10442
% de = f(e) + c + B(e)*v
% de ≈ A(e)e + B(e)v + (f(0) + c) 

% f(e) ≈ f(0) + A(e)e,
% Jacobian A from Eq. 16

n = 4;

% From [1]
% K = 1
% Q = 1000*eye(n)
% R = eye(n)
% dt = 0.01;
K = 1;
Q = 1000*eye(n);
R = eye(n);
dt = 0.01;

B = eye(n);
J = rand(n);
Ks = 0.1;

T = 0:dt:2;

% Reproduce Figure 2
% omega = [1.30, 1.39, 0.44, 1.28]';
% x_des = [-0.74, 0.27, 0.15]';
x0 = [0.60, 0.86, 0.84, -0.13]';


% x_des = rand(n,1);
x_des = randi([0 1],[n 1]);
% x = randi([0 1],[n 1]);

X = zeros(n,length(T));
U = zeros(n,length(T));
E = zeros(length(T), 1);
err = zeros(length(T), 1);
x = x0;
for t = 1:length(T)
    A = kuramotoSDC(x, K, Ks, J);

    [~, C] = icare(A, B, Q, R);

    u = -C*(x-x_des);
    
    dxdt = phaseModel(x, K, Ks, J);
    x = x + (dxdt+B*u)*dt;

    E(t) = lyapunov(x, K, Ks, J);
    err(t) = norm(x-x_des);

    X(:,t) = x;
    U(:,t) = u;
end


tiledlayout(3,1)
nexttile
plot(T, X)
ylabel('Phases (\pi)')
grid on
nexttile
plot(T, U)
xlabel('Time (s)');
ylabel('Control Input');
grid on;
nexttile
yyaxis left 
plot(T, E, '-x')
ylabel('Lyapunov')
hold on
yyaxis right
plot(T, err, '-r')
ylabel('Error')
xlabel('Time (cycles)');
grid on

function A = kuramotoSDC(x, K, Ks, J)
% Jacobian from state-dependent coefficients (SDC)
n = length(x);
A = zeros(n);
for i = 1:n
    for j = 1:n
        if i ~= j
            A(i,j) = K*pi*J(i,j)*cos(pi*(x(i)-x(j)));
        end
    end
    A(i,i) = -K*pi*sum(J(i,:).*cos(pi*(x(i)-x'))) - 2*Ks*cos(2*pi*x(i));
end
end