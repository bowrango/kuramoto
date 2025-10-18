function [X, Xrev, Xirrev, Xcheck] = isotropicOU2D(x0, L)
% Exact 2-dimensional multivariate Ornstein-Uhlenbeck process with
% isotropic diffusion. See Appendix A.2

d = 2;

theta = 0;
tau = 0.01;
eps = 0.5;

B = [2, -theta; theta, 2];
D = eps*eye(d);
S = [eps/2, 0; 0, eps/2];
Sinv = [2/eps, 0; 0, 2/eps];
expb = expm(-tau * B);
Q = B*S - D;
C = [exp(2*tau)*eps*sinh(2*tau), 0
    0, exp(2*tau)*eps*sinh(2*tau)];
expirrev = expm(-tau * Q * Sinv);
exprev = expm(-tau * D * Sinv);

X = zeros(d, L);
if nargout > 1
    Xrev = zeros(d, L);
    Xirrev = zeros(d, L);
    Xcheck = zeros(d, L);
    Xrev(:,1) = x0;
    Xirrev(:,1) = x0;
    Xcheck(:,1) = x0;
end

X(:,1) = x0;
for ii = 2:L
    mu = mvnrnd(zeros(d,1), C)';
    X(:,ii) = (expb * X(:,ii-1)) + mu;
    if nargout > 1
        Xrev(:,ii) = (exprev * Xrev(:,ii-1)) + mu;
        Xirrev(:,ii) = expirrev * Xirrev(:,ii-1);
        Xcheck(:,ii) = expirrev * (exprev * Xcheck(:,ii-1) + mu);
    end
end
end