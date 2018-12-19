% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, connFile] = ...
    connectEM.Consistency.loadConnectome(param);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Build axon-spine interface areas
asiT = ...
    connectEM.Connectome.buildAxonSpineInterfaces( ...
        param, graph, shAgglos, conn, syn);
asiT = asiT(asiT.type == 'PrimarySpine', :);
asiT.relId = reshape(1:height(asiT), [], 1);

%% Build connectivity matrix
clear cur*;

[pairT, ~, curIds] = unique(asiT(:, { ...
    'preAggloId', 'postAggloId'}), 'rows');
pairT.axonClass = conn.axonMeta.axonClass(pairT.preAggloId);
pairT.targetClass = conn.denMeta.targetClass(pairT.postAggloId);
pairT.asiIds = accumarray(curIds, asiT.relId, [], @(ids) {ids});

%% Actually analyse weight matrix
clear cur*;

%{
curAxonClasses = {'Corticocortical'};
curTargetClasses = {'OtherDendrite'};

curPairT = pairT( ...
    (isempty(curAxonClasses) | ismember(pairT.axonClass, curAxonClasses)) ...
  & (isempty(curTargetClasses) | ismember(pairT.targetClass, curTargetClasses)), :);

curPairT.medLog10AsiArea = cellfun( ...
    @(ids) median(log10(asiT.area(ids))), curPairT.asiIds);
%}

% Generate fake data
curDendCount = 10;
curAxonCount = 20;
curAxonSynCount = 100;

curAxonClasses = [ ...
    +0.5, -0.5, 1 / 3; ...
     0.0,  0.0, 1 / 3; ...
    -0.5, +0.5, 1 / 3];
curDendClasses = [1 / 2; 1 / 2];

rng(0);
curConn = cell(curAxonCount, curDendCount);

curAxonClassIds = repelem( ...
    1:size(curAxonClasses, 1), ceil(curAxonClasses(:, end)' * curAxonCount));
curAxonClassIds = curAxonClassIds(randperm(curAxonCount));
curDendClassIds = repelem( ...
    1:numel(curDendClasses), ceil(curDendClasses * curDendCount));
curDendClassIds = curDendClassIds(randperm(curDendCount));

for curAxId = 1:curAxonCount
    curAxClass = curAxonClassIds(curAxId);
    curDendIds = randi(curDendCount, 1, curAxonSynCount);
    
    curSynSizes = curDendClassIds(curDendIds);
    curSynSizes = curAxonClasses(curAxClass, curSynSizes);
    curSynSizes = curSynSizes + randn(size(curSynSizes));
    
    curConn(curAxId, :) = accumarray( ...
        curDendIds(:), curSynSizes(:), ...
       [curDendCount, 1], @(sizes) {sizes}, {zeros(0, 1)});
end

curMask = ~cellfun(@isempty, curConn);
curPairT = nan(sum(curMask(:)), 3);
[curPairT(:, 1), curPairT(:, 2)] = find(curMask);
curPairT(:, 3) = cellfun(@median, curConn(curMask));
curPairT = array2table(curPairT, 'VariableNames', { ...
    'preAggloId', 'postAggloId', 'medLog10AsiArea'});

[curAxonIds, ~, curPairT.preAggloId] = unique(curPairT.preAggloId);
[curDendIds, ~, curPairT.postAggloId] = unique(curPairT.postAggloId);

curAreas = [ ...
    curPairT.preAggloId, ...
    curPairT.postAggloId];
curAreas = accumarray( ...
    curAreas, curPairT.medLog10AsiArea, ...
   [numel(curAxonIds), numel(curDendIds)], [], nan);

Util.log('Clustering weight matrix');
[curLink, curDist] = connectEM.Consistency.linkage(curAreas);
[curDendLink, ~] = connectEM.Consistency.linkage(curAreas');
Util.log('Done!');
   
figure; [~, ~, curAxonPerm] = dendrogram(curLink, 0);
figure; [~, ~, curDendPerm] = dendrogram(curDendLink, 0);

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
curTemp = curDist(curAxonPerm, curAxonPerm);

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

xlim(curAx, [1, numel(curAxonIds)]);
xlabel(curAx, 'Axons (clustered)');
xticks(curAx, [1, numel(curAxonIds)]);
xticklabels(curAx, arrayfun( ...
    @num2str, xticks(curAx), 'UniformOutput', false));

ylim(curAx, [1, numel(curAxonIds)]);
ylabel(curAx, 'Axons (clustered)');
yticks(curAx, [1, numel(curAxonIds)]);
yticklabels(curAx, arrayfun( ...
    @num2str, yticks(curAx), 'UniformOutput', false));

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; 'Distance matrix'}, ...
    'FontWeight', 'normal', 'FontSize', 10);
