function plotRate(T, tess)
numStates = size(T, 1);
centroids = incenter(tess);
neighbors = tess.neighbors;

triplot(tess, 'k');
hold on

% Arrows: [x y dx dy intensity]
V = [];
for i = 1:numStates
    ci = centroids(i, :);
    for k = 1:3
        j = neighbors(i, k);
        if j > 0
            qji = T(j, i);
            if qji > 0
                cj = centroids(j, :);
                dir = cj - ci;
                dir = dir / norm(dir);
                V(end+1,:) = [ci dir qji];
            end
        end
    end
end

vecT = T(:);
xMax = max(vecT);
climit = [0, xMax];

cmap = parula;
numC = 256;
scale = 0.3;
for k = 1:size(V, 1)
    cidx = round(1 + (numC-1)*(V(k,5)-climit(1)) / (climit(2)-climit(1)));
    cidx = max(min(cidx, numC),1); % clamp
    color = cmap(cidx, :);
    quiver(V(k,1), V(k,2), V(k,3), V(k,4), scale, Color=color, MaxHeadSize=1, LineWidth=1);
end

colormap(cmap);
clim(climit);
cb = colorbar;
ylabel(cb, 'Intensity: $$q_{ji}$$', Interpreter='latex', FontSize=12);
title('Transition Vector Field')
% axis equal
xlabel('x1')
ylabel('x2')
% xlim([xMin xMax])
% ylim([xMin xMax])
end