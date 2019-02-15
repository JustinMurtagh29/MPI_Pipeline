% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = fullfile(curDir, sprintf('%s_asiT.mat', curAsiFile));

asiT = load(curAsiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Build connectivity matrix
clear cur*;
[~, curIds, curUniIds] = unique( ...
    asiT(:, {'preAggloId', 'postAggloId'}), 'rows');

pairT = asiT(curIds, ...
    {'preAggloId', 'postAggloId', 'axonClass', 'targetClass'});
pairT.asiIds = accumarray( ...
    curUniIds, reshape(1:height(asiT), [], 1), [], @(ids) {ids});

%% Actually analyse weight matrix
clear cur*;

curAxonClasses = {'Corticocortical'};
curTargetClasses = {'OtherDendrite'};

curPairT = pairT( ...
    (isempty(curAxonClasses) | ismember(pairT.axonClass, curAxonClasses)) ...
  & (isempty(curTargetClasses) | ismember(pairT.targetClass, curTargetClasses)), :);

[curAxonIds, ~, curPairT.preAggloId] = unique(curPairT.preAggloId);
[curDendIds, ~, curPairT.postAggloId] = unique(curPairT.postAggloId);

curPairT.medLog10AsiArea = cellfun( ...
    @(ids) median(log10(asiT.area(ids))), curPairT.asiIds);

areaMat = [ ...
    curPairT.preAggloId, ...
    curPairT.postAggloId];
areaMat = accumarray( ...
    areaMat, curPairT.medLog10AsiArea, ...
   [numel(curAxonIds), numel(curDendIds)], [], nan);

%% Actually do the clustering
Util.log('Clustering weight matrix');
[axonLink, axonDist] = connectEM.Consistency.linkage(areaMat);
[dendLink, dendDist] = connectEM.Consistency.linkage(areaMat');
Util.log('Done!');

curFig = figure;
[~, ~, axonPerm] = dendrogram(axonLink, 0);
[~, ~, dendPerm] = dendrogram(dendLink, 0);
close(curFig);
clear curFig;

%% Plot synapse sizes
clear cur*;

curPlots(1).data = areaMat;
curPlots(1).axonLink = [];
curPlots(1).dendLink = [];

curPlots(2).data = areaMat(axonPerm, dendPerm);
curPlots(2).axonLink = axonLink;
curPlots(2).dendLink = dendLink;

for curPlot = curPlots
    curData = curPlot.data;
    
    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [600, 600];
    
    curAx = axes(curFig); %#ok
    curIm = imagesc(curAx, curData);
    curIm.AlphaData = 1 - isnan(curData);
    curAx.CLim = [-1.5, 0.5];
    curAx.Color = 'black';
    curAx.Box = 'off';

    curCbar = colorbar('peer', curAx);
    curCbar.TickDirection = 'out';

    axis(curAx, 'equal');
    xlim(curAx, [1, size(curData, 2)]);
    ylim(curAx, [1, size(curData, 1)]);

    curAx.YDir = 'normal';
    curAx.Box = 'off';
    xticks(curAx, []);
    yticks(curAx, []);

    xlabel(curAx, 'Dendrites');
    xticks(curAx, curAx.XLim);
    xticklabels(curAx, arrayfun( ...
        @num2str, xticks(curAx), 'UniformOutput', false));
    curAx.XAxisLocation = 'top';

    ylabel(curAx, 'Axons');
    yticks(curAx, curAx.YLim);
    yticklabels(curAx, arrayfun( ...
        @num2str, yticks(curAx), 'UniformOutput', false));
    curAx.YAxisLocation = 'right';

    title(curAx, ...
        {info.filename; info.git_repos{1}.hash; 'Weight matrix'}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    curAx.Position = [0.15, 0.15, 0.65, 0.65];

    % Dendrite dendrogram
    if ~isempty(curPlot.dendLink)
        curTempFig = figure; 
        dendrogram(curPlot.dendLink, 0, 'Orientation', 'bottom');
        curTempFig.Children.XAxis.Visible = 'off';
        curTempFig.Children.YAxis.Visible = 'off';
        set(curTempFig.Children.Children, 'Color', 'black');

        copyobj(curTempFig.Children, curFig);
        close(curTempFig);

        curFig.Children(1).Position = [ ...
            curFig.Children(end).Position(1), 0, ...
            curFig.Children(end).Position([3, 2])];
    end

    % Axon dendrogram
    if ~isempty(curPlot.axonLink)
        curTempFig = figure; 
        dendrogram(curPlot.axonLink, 0, 'Orientation', 'left');
        curTempFig.Children.XAxis.Visible = 'off';
        curTempFig.Children.YAxis.Visible = 'off';
        set(curTempFig.Children.Children, 'Color', 'black');

        copyobj(curTempFig.Children, curFig);
        close(curTempFig);

        curFig.Children(1).Position = [ ...
            0, curFig.Children(end).Position(2)...
            curFig.Children(end).Position([1, 4])];
    end
end
