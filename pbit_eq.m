
% Stationary p-bit distribution

% [1] https://arxiv.org/pdf/2302.06457
% [2] https://arxiv.org/abs/2007.12331
% [3] https://arxiv.org/pdf/2503.10302
% [4] https://arxiv.org/pdf/2507.07420

% See sota [3];

% TODO What can I say for positive definite Q?

n = 8;
N = 2^n;
tau0 = 5;
tau = @(t) tau0/t;

rng default

% symmetric i.i.d Gaussian Qij ~ N(0, 1)
Q = randn(n);
h = diag(Q);
Q = tril(Q,-1) + triu(Q',1);

% 1. Markov Transition
% P1 = unitrans(Q, h, beta);
% mc = dtmc(P1);
% [probm, tMix] = asymptotics(mc);
% probm = probm';   
% [ev, lambda] = eigs(P1.', 2);
% sgap = 1 - lambda(2);
% epsilon = 1e-3;
% tmix = log(1/epsilon) / sgap;


T = 20;

% Boltzmann
% E(k) = -(0.5*X(:,k)'*Q*X(:,k) + h'*X(:,k));
X = (2*(dec2bin(0:N-1, n)=='1') - 1)';
Eb = -(sum(X.*(Q*X), 1)/2 + h'*X);
[Eb, perm] = sort(Eb);
% expE = exp(-beta.*Eb);
% probb = expE/sum(expE);
% probb = probb(perm);

numGibbsSteps = 1e3;

pbits = sign(randn(n, 1));
Et = zeros(numGibbsSteps,T);
F = zeros(N, T);
% counts = zeros(N, 1);
for t = 1:T
    beta = 1/tau(t);

    % Markov Transition
    Pt = unitrans(Q, h, beta);
    mc = dtmc(Pt);
    [probm, tMix] = asymptotics(mc);
    probm = probm';
    probm = probm(perm);
    
    % Gibbs Sampling
    counts = zeros(N, 1);
    for k = 1:numGibbsSteps
        % for ii = 1:n % system-scan
        %     I = W(ii,:)*pbits + h(ii);
        %     pbits(ii) = sign(tanh(beta*I) - (2*rand-1));
        % end
        idx = randi(n); % random-scan
    
        I = Q(idx,:)*pbits + h(idx);
        pbits(idx) = sign(tanh(beta*I) - (2*rand-1));
    
        bits = (pbits+1)/2;
        idx = bin2dec(join(string(bits'),""))+1;
        counts(idx) = counts(idx)+1;
        Et(k,t) = -(0.5*pbits'*Q*pbits + h'*pbits);
    end
    probg = counts/sum(counts);

    F(:,t) = probg(perm);
end

% aggregate probability by energy
[Ebin, ~, eidx] = unique(Eb);
Fbin = zeros(length(Ebin), T);
for t = 1:T
    for k = 1:length(Ebin)
        Fbin(k,t) = sum(F(eidx==k, t));
    end
end

% remove near-zero probability
minProb = 1/numGibbsSteps;
mask = any(Fbin >= minProb, 2);
Eplot = Ebin(mask);
Fplot = Fbin(mask,:);

[Tgrid, Egrid] = meshgrid(1:T, Eplot);
mask = Fplot >= minProb;
scatter(Tgrid(mask), Egrid(mask), 200, Fplot(mask), 'filled');
colormap(sky)
cbar = colorbar; 
ylabel(cbar,'Probability')
hold on

% TODO Misleading because +- std can go below the min(Eb)
% Emean = mean(Et);
% Estd = std(Et);
% errorbar(1:T, Emean, Estd, 'k.', 'LineWidth',1.5, 'MarkerSize',10);

yline(min(Eb))

xlabel('Time')
ylabel('Energy')
grid on


% hdl = imagesc(1:T, Eplot, Fplot);
% xticks(1:T)
% set(gca,'YDir','normal')
% colorbar
% xlabel('t')
% ylabel('Energy')

% hdl = bar3(Eplot, Fplot, 'blue');
% xlabel('t')
% ylabel('Energy')
% zlabel('Probability')

% for k = 1:numel(hdl)
%     Z = hdl(k).ZData;
%     mask = Z < minProb;
%     hdl(k).AlphaData = ~mask;
%     hdl(k).FaceAlpha = 'flat';
%     hdl(k).EdgeColor = 'none';
% end

% Only for large gibbs steps
% TV = 0.5 * sum(abs(probg - probmc));
% assert(TV < epsilon);

% Kmax = 50;                 
% tv = zeros(Kmax,1);
% % worst-case initial state
% mu = zeros(N,1);
% mu(1) = 1;
% % tv(k) decays like exp(-k/tmix)
% for k = 1:Kmax
%     tv(k) = 0.5 * sum(abs(mu - probm));
%     mu = T1.'*mu;
% end
% figure;
% semilogy(1:Kmax, tv, 'o-');
% grid on;
% xlabel('Single-spin updates (steps of T1)');
% ylabel('Total variation distance');

% [Eb, perm] = sort(Eb);
% probg = probg(perm);
% probb = probb(perm);
% probm = probm(perm);

% figurex
% spy(T)
% title('P(x->y)')
% label = 1:N;
% xticks(label)
% xticklabels(dec2bin(label-1))
% yticks(label)
% yticklabels(dec2bin(label-1))

% TODO yyaxis quality
% figure
% grid on
% hold on
% plot(Eb, probg, 'bx')
% plot(Eb, probb, 'ro')
% plot(Eb, probm, 'r-')
% xlabel('Energy H(x)')
% ylabel('Probability')
% legend('Gibbs', 'Boltzmann', 'Markov', Location='best')


% Compare to SK theory

% sigma2 = sum(triu(Q,1).^2, 'all') + sum(h.^2);
% 
% betas = linspace(0, 1.0, 40);
% E_emp = zeros(size(betas));     
% E_the = zeros(size(betas));    
% for ii = 1:length(betas)
%     beta = betas(ii);
% 
%     w = exp(-beta*E);
%     w = w/sum(w);
% 
%     E_emp(ii) = sum(w .* E);
%     E_the(ii) = -beta * sigma2;
% end
% 
% figure; hold on; grid on;
% plot(betas, E_emp/n, 'bo-', 'LineWidth', 1.2, 'DisplayName', 'Empirical');
% plot(betas, E_the/n, 'r--', 'LineWidth', 1.5, 'DisplayName', 'SK Theory');
% 
% xlabel('\beta');
% ylabel('Average Spin Energy <E>/n');
% legend('Location','best');


function P = unitrans(Q, h, beta)
% Transition matrix for a uniform random p-bit update (random scan)
% P(i,j) is the probability state i transitions to state j.
% Rows sum to 1. States with Hamming distance more than 1 have no
% probability.
n = size(Q,1);
N = 2^n;
X = 2*(dec2bin(0:N-1, n)=='1') - 1;
P = zeros(N);
for ii = 1:N
    x = X(ii, :);
    for idx = 1:n
        jj = bitxor(ii-1, bitshift(1, n-idx)) + 1;

        f = (Q(idx,:)*x' + h(idx));
        % p = (1 - x(idx)*tanh(beta*f))/2;    
        p = 1 / (1 + exp(2*beta*x(idx)*f));
        
        P(ii,jj) = P(ii,jj) + p/n;
        P(ii,ii) = P(ii,ii) + (1-p)/n;
    end
end
end

function P = systrans(W, h, beta)
% Transition matrix for an all p-bit update (system scan)
% P(i,j) is the probability state i transitions to state j.
% Rows sum to 1. All elements are non-zero.
% FIXME negative eigs
n = size(W, 1);
N = 2^n;
states = 2*(dec2bin(0:N-1, n)=='1') - 1;
P = zeros(N);
for i = 1:N
    s = states(i, :)';
    f = W*s + h;
    for j = 1:N
        sp = states(j, :)';
        p = 0.5 * (1 + sp .* tanh(beta*f));
        P(i,j) = prod(p);
    end
end
end