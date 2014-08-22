function visualizeEdgeFeaturesNewest(pT)

load(pT.gp.initalGroundTruth);
% Create output directory
if ~exist('/zdata/manuel/sync/feature/');
    mkdir('/zdata/manuel/sync/feature/');
end
% Visualizierungs-Settings
colors = {'g' 'r' 'g' 'r'};
lineStyle = {'-' '-' ':' ':'};
% Training Daten aufsplitten
isEdgeTR = trainingData(trainingLabels == 1,:);
isNonEdgeTR = trainingData(trainingLabels == -1,:);
isEdgeTE = testData(testLabels == 1,:);
isNonEdgeTE = testData(testLabels == -1,:);
% Titel vorberechnen
titleString = [num2str(size(isEdgeTR,1)) ' training positive, ' num2str(size(isNonEdgeTR,1)) ' training negative, ' ...
    num2str(size(isEdgeTE,1)) ' test positive, ' num2str(size(isNonEdgeTE,1)) ' test negative, ' ];

for idx=1:size(isEdgeTR,2)
    % Bin Vektor berrechnen
    limits = [min([isEdgeTR(:,idx); isNonEdgeTR(:,idx); isEdgeTE(:,idx); isNonEdgeTE(:,idx)]) ...
        max([isEdgeTR(:,idx); isNonEdgeTR(:,idx); isEdgeTE(:,idx); isNonEdgeTE(:,idx)])];
    binSize = (limits(2) - limits(1))/100;
    x = limits(1)+binSize/2:binSize:limits(2)-binSize/2;

    a = hist(isEdgeTR(:,idx),x);
    a = a ./ sum(a);
    b = hist(isNonEdgeTR(:,idx),x);
    b = b ./ sum(b);

    e = hist(isEdgeTE(:,idx),x);
    e = e ./ sum(e);
    f = hist(isNonEdgeTE(:,idx),x);
    f = f ./ sum(f);

    figure('Visible', 'off', 'Units', 'centimeters', 'Position', [0 0 29.7 21], ...
        'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
    h = stairs(x, [a; b; e; f]');
    for j=1:length(h)
        set(h(j),'Color',colors{j}, 'LineStyle', lineStyle{j}, 'LineWidth', 2);
    end
    title(titleString);
    legend('training connected', 'training unconnected', 'test connected', 'test unconnected');
    set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
    print(gcf, ['/zdata/manuel/sync//feature/' num2str(idx, '%.3i') '.pdf'], '-dpdf');
    close all;
end

end

