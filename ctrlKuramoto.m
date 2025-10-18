
close all

% There can be muliple stable non-global stationary distributions. How do 
% we guide the system with u,K,Ks? Kuramoto as SDE Decompose station distribution

rng(0)

% FIXME control gain
% TODO natural frequencies
% TODO B matrix

% [1] https://arxiv.org/pdf/2504.18399
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
omega = [1.30, 1.39, 0.44, 1.28]';
dx_des = [-0.74, 0.27, 0.15]';
x0 = [0.60, 0.86, 0.84, -0.13]';


% x_des = rand(n,1);
% x_des = randi([0 1],[n 1]);
% x = randi([0 1],[n 1]);

X = zeros(n,length(T));
U = zeros(n,length(T));
E = zeros(length(T), 1);
err = zeros(length(T), 1);
for t = 1:length(T)
    A = kuramotoSDC(x, K, Ks, J);

    [~, C] = icare(A, B, Q, R);

    u = -C*(x-x_des);
    
    dxdt = phaseModel(x, K, Ks, J);
    x = x + (dxdt+B*u)*dt;

    E(t) = energy(x, K, Ks, J);
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

function dxdt = phaseModel(x, K, Ks, J)
n = length(x);
dxdt = zeros(n,1);
for i = 1:n
    dxdt(i) = -K*(J(i,:)*sin(pi*(x(i)-x)));
end
dxdt = (dxdt - Ks*sin(2*pi*x))/pi;
end

function A = kuramotoSDC(x, K, Ks, J)
% Jacobian
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

function E = energy(x, K, Ks, J)
% Lyapunov (4.7)
shift = x - x.';
E = -K*sum(J(:).*cos(pi*(shift(:)))) - Ks*sum(cos(2*pi*x));
end