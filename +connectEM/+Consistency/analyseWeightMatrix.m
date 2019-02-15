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

curAreas = [ ...
    curPairT.preAggloId, ...
    curPairT.postAggloId];
curAreas = accumarray( ...
    curAreas, curPairT.medLog10AsiArea, ...
   [numel(curAxonIds), numel(curDendIds)], [], nan);

%% Actually do the clustering
Util.log('Clustering weight matrix');
[curAxonLink, curAxonDist] = connectEM.Consistency.linkage(curAreas);
[curDendLink, curDendDist] = connectEM.Consistency.linkage(curAreas');
Util.log('Done!');

curFig = figure;
[~, ~, curAxonPerm] = dendrogram(curAxonLink, 0);
[~, ~, curDendPerm] = dendrogram(curDendLink, 0);
close(curFig);
clear curFig;

%% Show dendrogram
curPlots = struct;
curPlots(1).type = 'Axons'; curPlots(1).link = curAxonLink;
curPlots(2).type = 'Dendrites'; curPlots(2).link = curDendLink;

for curPlot = curPlots
    curFig = figure();
    curFig.Color = 'white';
    
    curAx = axes(curFig); %#ok
    dendrogram(curPlot.link, 0);
    xlabel(curAx, curPlot.type);
    xticks(curAx, []);
    curAx.TickDir = 'out';
    
    title(curAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

%%
%{
curDataFile = '/mnt/mpibr/data/Personal/mottaa/L4/2018-12-19-Clustering-of-Weight-Matrix/corticocortical-onto-other-dendrites.mat';
curData = load(curDataFile);
curAxonLink = curData.curAxonLink;
curAxonLink = curAxonLink(not(isnan(curAxonLink(:, end))), :);
%}

%% Find threshold with a second reasonably-sized cluster
assert(issorted(curAxonLink(:, 1:2), 2));

curAxonCount = max(reshape(curAxonLink(:, 1:2), [], 1));
curAxonCount = curAxonCount + 1 - size(curAxonLink, 1);

curOut = nan(0, 3);

for i = 1:500
    curGraph = graph( ...
        cat(1, curAxonLink(:, 1), curAxonLink(:, 2)), ...
        repmat(transpose((1:size(curAxonLink, 1)) + curAxonCount), 2, 1));

    curCompIds = conncomp(curGraph);
    curCompSizes = accumarray(curCompIds(:), 1);
    
    curSecCompSize = feval( ...
        @(s) reshape(s(1:2), 1, []), ...
        sort([0; curCompSizes], 'descend'));
    curOut(end + 1, :) = [size(curOut, 1), curSecCompSize]; %#ok
    
    curAxonLink(end, :) = [];
end

% Plot
curFig = figure;
curFig.Color = 'white';
curFig.Position(3:4) = [310, 180];
curAx = axes(curFig);
hold(curAx, 'on');

plot(curOut(:, 1), curOut(:, 2), 'LineWidth', 2);
plot(curOut(:, 1), curOut(:, 3), 'LineWidth', 2);

set(curAx, 'TickDir', 'out', 'YScale', 'log');

curLeg = legend(curAx, { ...
    'Largest cluster', 'Second largest cluster'});
set(curLeg, 'Location', 'SouthEast', 'Box', 'off');

xlabel(curAx, 'Linkage cut-off');
ylabel(curAx, 'Cluster size');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Extract the largest two clusters
curAxonLink = curData.curAxonLink;
curAxonLink = curAxonLink(not(isnan(curAxonLink(:, end))), :);
assert(issorted(curAxonLink(:, 1:2), 2));

curAxonCount = max(reshape(curAxonLink(:, 1:2), [], 1));
curAxonCount = curAxonCount + 1 - size(curAxonLink, 1);

% Use linkage cut-off of 200
curAxonLink((end - 199):end, :) = [];

curGraph = graph( ...
    cat(1, curAxonLink(:, 1), curAxonLink(:, 2)), ...
    repmat(transpose((1:size(curAxonLink, 1)) + curAxonCount), 2, 1));

curCompIds = conncomp(curGraph);
curCompSizes = accumarray(curCompIds(:), 1);
[~, curSortIds] = sort(curCompSizes, 'descend');

curAxonIdsA = find(curCompIds(1:curAxonCount) == curSortIds(1));
curAxonIdsB = find(curCompIds(1:curAxonCount) == curSortIds(2));

%% Plot synapse sizes
curTemp = curAreas(curAxonPerm, curDendPerm);

curFig = figure();
curFig.Color = 'white';
curAx = axes(curFig);
curIm = imagesc(curAx, curTemp);
curIm.AlphaData = 1 - isnan(curTemp);
curAx.Color = 'black';
curAx.Box = 'off';

curCbar = colorbar('peer', curAx);
curCbar.TickDirection = 'out';
curCbar.Label.String = 'log_{10}(median axon-spine interface area [µm²])';

axis(curAx, 'equal');
xlim(curAx, [1, numel(curDendIds)]);
ylim(curAx, [1, numel(curAxonIds)]);

xlabel(curAx, 'Dendrites');
ylabel(curAx, 'Axons (clustered)');

curAx.TickDir = 'out';
xticks(curAx, [1, numel(curDendIds)]);
xticklabels(curAx, arrayfun( ...
    @num2str, xticks(curAx), 'UniformOutput', false));
yticks(curAx, [1, numel(curAxonIds)]);
yticklabels(curAx, arrayfun( ...
    @num2str, yticks(curAx), 'UniformOutput', false));

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; 'Weight matrix'}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot distance matrix
curIds = [curAxonIdsA, curAxonIdsB];
curTemp = curAxonDist(curIds, curIds);

curFig = figure();
curFig.Color = 'white';
curAx = axes(curFig);
% NOTE(amotta): For some reason MATLAB throws an error if `curAx` is passed
% as first input argument to `imagesc` (despite this being a valid input
% according to the documentation).
curIm = imagesc(curTemp, [0, 1.14]); % HACK(amotta): Hard-coded limits!
curIm.AlphaData = 1 - isnan(curTemp);
axis(curAx, 'square');
curAx.Color = 'black';
curAx.Box = 'off';

curCbar = colorbar('peer', curAx);
curCbar.TickDirection = 'out';
curCbar.Label.String = 'Squared Euclidean distance';

axis(curAx, 'square');
curAx.TickDir = 'out';
curAx.XAxisLocation = 'top';

% xlim(curAx, [1, numel(curAxonIds)]);
xlabel(curAx, 'Axons (clustered)');
% xticks(curAx, [1, numel(curAxonIds)]);
xticklabels(curAx, arrayfun( ...
    @num2str, xticks(curAx), 'UniformOutput', false));

% ylim(curAx, [1, numel(curAxonIds)]);
ylabel(curAx, 'Axons (clustered)');
% yticks(curAx, [1, numel(curAxonIds)]);
yticklabels(curAx, arrayfun( ...
    @num2str, yticks(curAx), 'UniformOutput', false));

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; 'Distance matrix'}, ...
    'FontWeight', 'normal', 'FontSize', 10);
