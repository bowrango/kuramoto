
nOsc = 100;

rng default

% Random problem graph.
% (See https://arxiv.org/pdf/1811.11538)
% resolution = 2;

A = abs(randn(nOsc));
A = (A+A.')/2;
A(1:nOsc+1:end) = 0;

% A = [ 0 0.9294 0.1682 0.2574 0.1395 0.3284
%       0.9294 0 0.8639 0.8428 1.0879 0.0015
%       0.1682 0.8639 0 1.2203 1.3125 0.0177
%       0.2574 0.8428 1.2203 0 0.8108 0.5174
%       0.1395 1.0879 1.3125 0.8108 0 1.6862
%       0.3284 0.0015 0.0177 0.5174 1.6862 0];

% A = randi([1 nOsc], nOsc);
% A = (A+A.')/2;
% A(1:nOsc+1:end) = 0;

% A = [0     2     0     0   -2.2   0     0     3   ;
%     2     0     0.8   0    0     2     0     0   ;
%     0     0.8   0     1.3  0     0     1.2   0   ;
%     0     0     1.3   0    2     0     0     0.2 ;
%    -2.2   0     0     2    0    -0.5   0     0   ;
%     0     2     0     0   -0.5   0     0.5   0   ;
%     0     0     1.2   0    0     0.5   0     1   ;
%     3     0     0     0.2  0     0     1     0   ];

% A = [0 0 1 0;
%      0 0 1 0;
%      1 1 0 1;
%      0 0 1 0];

% A = randi(2^resolution, nOsc);
% A(1:nOsc+1:end) = 0;
% A = (A+A.');
% 
% A = [0 1 2;
%      1 0 1
%      2 1 0];

% A = ones(nOsc);
% A(1:nOsc+1:end) = 0;

% A = [0 1 0 0
%      1 0 1 0
%      0 1 0 1
%      0 0 1 0];
g = graph(A);

% Find groundtruth with classical solver.
% max sum{aij*(xi + xj - 2*xi*xj)} -> min x.'*Q*x + r.'*x
qb = maxcut2qubo(g);
% qb = qubo(-A);
gt = solve(qb);

% Compute MaxCut. Consider nOsc=2:1
% When Aij > 0, bits are opposite because we cut the only edge
% When Aij < 0, bits are 0 because no cut is made. The MaxCux is 0.
bestx = gt.BestX;
maxcut = (~bestx.')*A*bestx;

% Formulate MaxCut as Ising
J = -A;

% L = J ./ abs(sum(J, "all"));

% TODO device frequency
tstop = 2;
dt = 1e-2;

% TODO K contribute to number of cycles

% Coupling
% K = 1;
K = 3;
carg.k = K/tstop;
carg.T = tstop/K;
carg.K = K;
coup = @(t, args) t*args.k;
% coup = @(t, K) K*tanh(cos(pi*t));
% coup = @(t, arg) K.*t.*sin(t);
% coup = @(t, args) args.K*ones(size(t));
% coup = @(t, args) 7*tanh(cos(2*pi.*t./args.T));

% Sync
Ks = 5.3;
sarg.T = tstop/10;
sarg.Ks = Ks;
sarg.k = Ks/tstop;
sync = @(t, args) args.Ks*ones(size(t)); % const.
% sync = @(t, args) t.*args.k; % ramp
%sync = @(t, args)4+6*tanh(10*cos(2*pi.*t./args.T));

drift = @(t,X) phaseModel(X, coup(t,carg), sync(t,sarg), J);

% Noise
% TODO relation inverse temp.
Kn = 0.1;
diffusion = @(t,X) Kn*eye(nOsc);

% Initial phase (radians)
% rng(0)
x0 = rand(nOsc,1);

model = sde(drift, diffusion, StartState=x0);
[X, T] = simulate(model, round(tstop/dt), DeltaTime=dt);

E = zeros(size(T));
cuts = zeros(size(T));
ising = zeros(size(T));
order = zeros(size(T));
for k = 1:length(T)
    x = X(k,:)';
    E(k) = lyapunov(x, coup(T(k),carg), sync(T(k),sarg), J);
    % Ising function
    % Convert phases to spins
    s = ones(nOsc,1);
    odds = find(mod(round(x), 2));
    s(odds) = -1;
    % Convert to binary
    bits = (s+1)/2;
    % Compute MaxCut. Using spins its (s.'*J*s)/2
    cuts(k) = (~bits.'*A*bits);
    ising(k) = -(s.'*J*s)/2;
    % Order parameter
    % order(k) = abs(1/nOsc*sum(exp(1j*x)));
end

hdl = tiledlayout(4,1);
title(hdl, "n="+string(nOsc)+" oscillators, randn J_{ij}<0")

nexttile
plot(T, X)
xlim([0 tstop])
% title('No SYNC')
% ylabel('phases (\pi)')
ylabel('Phase (\pi)')
% legend("osc"+string(1:min([nOsc 8])))
grid on

nexttile
hold on
% Negate to make cuts positive
yline(-gt.BestFunctionValue, LineWidth=2)
plot(T, cuts, '-x')
%plot(T, ising, 'o')
xlim([0 tstop])
% ylabel('Ising')
ylabel('Cut')
legend('groundstate', 'value', Location='best')
grid on
hold off

nexttile
hold on
plot(T, E, '-x')
ylabel('Lyapunov')
xlim([0 tstop])
grid on
hold off

nexttile
hold on
grid on
yyaxis right
plot(T, coup(T,carg), '--')
xlim([0 tstop])
ylabel('K')
yyaxis left
plot(T, sync(T,sarg), '--')
xlim([0 tstop])
xlabel('Time (cycles)');
ylabel('K_{s}')
yyaxis left
hold off

% figure
% plot(T, order)

dir = "results/";
filename = replace(string(datetime('now')), " ", "-")+".png";
saveas(gcf, dir+filename)