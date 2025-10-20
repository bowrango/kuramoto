
% https://arxiv.org/pdf/2409.07479

rng default

d = 2;
numSamples = 100;
numPaths = 5;
numDownSamples = 100;

K = 1;
coup = @(t, K) K;
Ks = 5;
sync = @(t, Ks) Ks;

% J = randn(d);
% J = (J+J.')/2;
J = eye(d);
drift = @(t,X) phaseModel(X, coup(t, K), sync(t, Ks), J);
Kn = 0.1;
diffusion = @(t,X) Kn*eye(d);

x0 = rand(d,1);
model = sde(drift, diffusion, StartState=x0);
X = simulate(model, numSamples, NTrials=numPaths);
X = permute(X, [2 1 3]);

% X = zeros(d,numSamples,numPaths);
% for ii=1:numPaths
%     X(:,:,ii) = isotropicOU2D(rand(1,d), numSamples);
% end

figure
plot(X(:,:,1)')
xlabel('step')
legend('x1', 'x2')
title('Feature Path')

% figure
% scatter(X(1,:,1), X(2,:,1))
% title('2A')

U = reshape(X, 2, []);
S = downsample(U, numDownSamples);
% figure
% scatter(S(1,:), S(2,:))
% title('2B')

DT1 = delaunayTriangulation(S(1,:)', S(2,:)');
IC1 = incenter(DT1);
figure
triplot(DT1)
hold on
plot(IC1(:,1),IC1(:,2),'*r')
title('2C')
hold off

% FIXME Generator matrix
Q = transitionRate(DT1, X);
% plotRate(Q, DT1)
% title('2D')

% figure
% subplot(2,1,1)
% plot(X(:,:,1)')
% xlabel('step')
% legend('x1', 'x2')
% title('Trajectory')
% subplot(2,1,2)
% plotRate(Q, DT1)



DT2 = delaunayTriangulation(IC1(:,1), IC1(:,2));
IC2 = incenter(DT2);
% figure
% triplot(DT2)
% hold on
% plot(IC2(:,1),IC2(:,2),'*r')
% title('2E')
% hold off

% Edge flow is analogous to the excess work required to move from
% state i to j (Eq. 10).
% TODO KL Divergence (rel entropy)
F = edgeflow(Q, DT2);

plotFlow(F, DT2, 'Flow')

G = gradmatrix(DT2);
y = G\F;
Fgrad = G*y;
Fcurl = F - Fgrad;

plotFlow(Fgrad, DT2, 'Grad')

plotFlow(Fcurl, DT2, 'Curl')

function G = gradmatrix(DT)
e = edges(DT);
G = incidence(graph(e(:,1), e(:,2)))';
end


function ef = edgeflow(Q, DT)
e = edges(DT);
ef = zeros(length(e),1);
for ii = 1:length(e)
    if Q(e(ii,1),e(ii,2)) > 0
        qij = Q(e(ii,1),e(ii,2));
        qji = Q(e(ii,2),e(ii,1));
        ef(ii) = log(sqrt(qij/qji));
    end
end
end

function [S, idx] = downsample(P, numPoints)
% Furthest point sampling implementation
% re: Brunton compressed sensing with random points
D = size(P,2);
idx = zeros(numPoints, 1);
distances = inf(1, D);

idx(1) = randi(D);  % Random starting point
for i = 2:numPoints
    diff = P - P(:,idx(i-1));
    dist = sqrt(sum(diff.^2, 1));
    distances = min(distances, dist);
    [~, idx(i)] = max(distances);
end
S = P(:, idx);
end