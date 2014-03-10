function visualizeEdgeFeaturesNewest( train, test, weightLabels, idx)
% Visualizierungs-Settings
colors = {'g' 'r' 'b' 'y'};
lineStyle = {':' ':' '-' '-'};
% Training Daten aufsplitten
isEdgeTR = train.X(train.y == 1,:);
isNonEdgeTR = train.X(train.y == -1,:);
% Bin Vektor berrechnen
limits = [min([isEdgeTR(:,idx); isNonEdgeTR(:,idx); test.X(:,idx)]) ...
    max([isEdgeTR(:,idx); isNonEdgeTR(:,idx); test.X(:,idx)])];
binSize = (limits(2) - limits(1))/40;
x = limits(1)+binSize/2:binSize:limits(2)-binSize/2;

a = hist(isEdgeTR(:,idx),x);
a = a ./ sum(a);
b = hist(isNonEdgeTR(:,idx),x);
b = b ./ sum(b);

e = hist(train.X(:,idx),x);
e = e ./ sum(e);
f = hist(test.X(:,idx),x);
f = f ./ sum(f);

figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
h = stairs(x, [a; b; e; f]');
for j=1:length(h)
    set(h(j),'Color',colors{j}, 'LineStyle', lineStyle{j}, 'LineWidth', 2);
end
title(weightLabels{idx});
legend('training connected', 'training unconnected', 'training set', 'test set');
set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
if ~exist('C:/Users/mberning/Desktop/classifier/features/');
    mkdir('C:/Users/mberning/Desktop/classifier/features/');
end
print(gcf, ['C:\Users\mberning\Desktop\classifier\features\feature' num2str(idx, '%.2i') '.pdf'], '-dpdf');
close all;
end
