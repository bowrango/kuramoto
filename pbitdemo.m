% Matt Bowring

% Compute the stationary distribution of p-bits using 3 methods.

% [1] https://arxiv.org/pdf/2302.06457
% [2] https://arxiv.org/abs/2007.12331
% [3] https://arxiv.org/pdf/2503.10302
% [4] https://arxiv.org/pdf/2507.07420

% See sota [3];

n = 10;
N = 2^n;
% beta = 0.1;
beta = 0.1;

rng default

W = randn(n);
W = (W + W')/2;
h = diag(W);
W(1:n+1:end)=0;

% Markov Transition
T1 = onetrans(W, h, beta);
T2 = systrans(W, h, beta);
mc = dtmc(T1);
[probmc,tMix] = asymptotics(mc);
probmc = probmc';   
[ev, lambda] = eigs(T1.', 2);
probev = ev(:,1) / sum(ev(:,1));
% gap = 1 - lambda(2);
% epsilon = 1e-3;
% tmix = log(1/epsilon) / gap;
% numGibbsSteps = ceil(n*tmix);
numGibbsSteps = 1e5;

% Gibbs sampling
pbits = sign(randn(n, 1));
counts = zeros(2^n, 1);
for t = 1:numGibbsSteps
    % for ii = 1:n % system-scan
    %     I = W(ii,:)*pbits + h(ii);
    %     pbits(ii) = sign(tanh(beta*I) - (2*rand-1));
    % end
    idx = randi(n); % random-scan
    I = W(idx,:) * pbits + h(idx);
    pbits(idx) = sign(tanh(beta * I) - (2*rand-1));
    bits = (pbits+1)/2;
    idx = bin2dec(join(string(bits'),""))+1;
    counts(idx) = counts(idx) + 1;
end
probg = counts/sum(counts);

% Only for large gibbs steps
% TV = 0.5 * sum(abs(probg - probmc));
% assert(TV < epsilon);

% Boltzmann
states = (2*(dec2bin(0:N-1, n)=='1') - 1)';
E = (sum(states.*(W*states), 1)/2 + h'*states)';
expE = exp(-beta.*E);
probb = expE/sum(expE);

% figure
% spy(T)
% title('P(x->y)')
% label = 1:N;
% xticks(label)
% xticklabels(dec2bin(label-1))
% yticks(label)
% yticklabels(dec2bin(label-1))

% TODO yyaxis quality
figure
grid on
hold on
plot(sort(probg), 'r.')
plot(sort(probb), 'm-')
% plot(sort(probev), 'g-')
plot(sort(probmc), 'k-')
xlim([0 N])
xlabel('State')
ylabel('Probability')
legend('Gibbs', 'Boltzmann', 'dtmc', Location='best')

function P = onetrans(W, h, beta)
% Transition matrix for one uniform random p-bit update (random scan)
% P(i,j) is the probability state i transitions to state j.
% Rows sum to 1. Zeros are on impossible states.
n = size(W,1);
N = 2^n;
states = 2*(dec2bin(0:N-1, n)=='1') - 1;
P = zeros(N);
for ii = 1:N
    m1 = states(ii, :);
    for idx = 1:n
        jj = bitxor(ii-1, bitshift(1, n-idx)) + 1;
        m2 = m1;
        m2(idx) = -m1(idx);
        f = W(idx,:)*m2' + h(idx);
        p = (m2(idx)*tanh(beta*f) + 1) / 2;
        P(ii,jj) = P(ii,jj) + p/n;
        P(ii,ii) = P(ii,ii) + (1-p)/n;
    end
end
end

function P = systrans(W, h, beta)
% Transition matrix for all p-bit updates (system scan)
% P(i,j) is the probability state i transitions to state j.
% Rows sum to 1. All elements are non-zero.
n = size(W, 1);
N = 2^n;
states = 2*(dec2bin(0:N-1, n)=='1') - 1;
P = zeros(N);
for i = 1:N
    s = states(i, :)';
    f = W*s + h;
    for j = 1:N
        sp = states(j, :)';
        p = 0.5 * (1 + sp .* tanh(beta * f));
        P(i,j) = prod(p);
    end
end
end