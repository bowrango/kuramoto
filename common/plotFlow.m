function plotFlow(ef, DT, label)

Y = DT.Points;
E = edges(DT);

cmap = parula;
nColors = size(cmap, 1);

minFlow = min(ef);
maxFlow = max(ef);

maxflow = max([maxFlow /2, -minFlow/2]);

figure;
triplot(DT, 'k');
hold on;
axis equal;

for i = 1:size(E, 1)
    e = E(i,:);
    if ef(i) < 0
        p1 = Y(e(1),:);
        p2 = Y(e(2),:);
    elseif ef(i) > 0
        p1 = Y(e(2),:);
        p2 = Y(e(1),:);
    else
        continue;
    end

    flowval = abs(ef(i));
    
    % Normalize flow to colormap index
    cidx = round((flowval / maxflow) * (nColors - 1)) + 1;
    cidx = min(max(cidx, 1), nColors);  % Clamp to valid index
    color = cmap(cidx, :);
    
    d = p2 - p1;  % direction vector
    quiver(p1(1), p1(2), d(1), d(2), Color=color, MaxHeadSize=0.05, LineWidth=1.2);
end

colormap(cmap);
clim([minFlow maxFlow]);
cb = colorbar;
ylabel(cb, 'Edge Flow', Interpreter='latex', FontSize=12);
if nargin < 3
    title('Edge Flow');
else
    title(label)
end