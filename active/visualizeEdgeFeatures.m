function visualizeEdgeFeatures( train, test, weightLabels, saveLocation)

colors = {'g' 'r' 'y'};
lineStyle = {'-' '-' ':'};
train.isEdge = train.X(:,train.y == 1);
train.isNonEdge = train.X(:,train.y == 0);

for i=1:length(weightLabels)
    figure('Visible', 'off');
    limits = [min([train.isEdge(i,:) train.isNonEdge(i,:) test.X(i,:)]) ...
        max([train.isEdge(i,:) train.isNonEdge(i,:) test.X(i,:)])];
    binSize = (limits(2) - limits(1))/50;
    x = limits(1)+binSize/2:binSize:limits(2)-binSize/2;
    a = hist(train.isEdge(i,:),x);
    a = a ./ sum(a);
    b = hist(train.isNonEdge(i,:),x);
    b = b ./ sum(b);
    c = hist(test.X(i,:),x);
    c = c ./ sum(c);
    h = stairs(x, [a; b; c]');
    for j=1:length(h)
        set(h(j),'Color',colors{j}, 'LineStyle', lineStyle{j}, 'LineWidth', 2);
    end
    title(weightLabels{i});
    legend('training set connected', 'training set unconnected', 'test set');
    set(gcf, 'PaperPositionMode', 'manual', 'PaperType', 'A4', 'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
    if ~exist([saveLocation 'features/'], 'dir')
        mkdir([saveLocation 'features/'])
    end
    saveas(gcf, [saveLocation 'features/' genvarname(weightLabels{i}) '.pdf']);
    close all;
end



end
