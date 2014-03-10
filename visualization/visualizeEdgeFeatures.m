function visualizeEdgeFeatures( trainSk, trainVol, test, weightLabels, idx)
% Visualizierungs-Settings
colors = {'g' 'r' 'g' 'r' 'b' 'y' 'c'};
lineStyle = {':' ':' '--' '--' '-' '-' '-'};
% Skeleton Daten aufsplitten
isEdgeS = trainSk.X(trainSk.y == 1,:);
isNonEdgeS = trainSk.X(trainSk.y == -1,:);
% Volume Daten aufsplitten
isEdgeV = trainVol.X(trainVol.y == 1,:);
isNonEdgeV = trainVol.X(trainVol.y == -1,:);
% Bin Vektor berrechnen
limits = [min([isEdgeS(:,idx); isNonEdgeS(:,idx); isEdgeV(:,idx); isNonEdgeV(:,idx); test(:,idx)]) ...
    max([isEdgeS(:,idx); isNonEdgeS(:,idx); isEdgeV(:,idx); isNonEdgeV(:,idx); test(:,idx)])];
binSize = (limits(2) - limits(1))/100;
x = limits(1)+binSize/2:binSize:limits(2)-binSize/2;

a = hist(isEdgeS(:,idx),x);
a = a ./ sum(a);
b = hist(isNonEdgeS(:,idx),x);
b = b ./ sum(b);

c = hist(isEdgeV(:,idx),x);
c = c ./ sum(c);
d = hist(isNonEdgeV(:,idx),x);
d = d ./ sum(d);

e = hist(trainSk.X(:,idx),x);
e = e ./ sum(e);
f = hist(trainVol.X(:,idx),x);
f = f ./ sum(f);
g = hist(test(:,idx),x);
g = g ./ sum(g);

figure('position', [1 41 1600 784], 'Renderer', 'OpenGL');
h = stairs(x, [a; b; c; d; e; f; g]');
for j=1:length(h)
    set(h(j),'Color',colors{j}, 'LineStyle', lineStyle{j}, 'LineWidth', 2);
end
title(weightLabels{idx});
legend('skeleton training connected', 'skeleton training unconnected', 'volume training connected', 'volume training unconnected', 'skeleton training set', 'volume training set', 'test set');
set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
saveas(gcf, ['C:/Users/mberning/Desktop/feature/featureNumber' num2str(idx, '%.3i') '.pdf']);
close all;

end
