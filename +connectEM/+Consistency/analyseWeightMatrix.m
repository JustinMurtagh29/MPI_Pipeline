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
asiT(asiT.type ~= 'PrimarySpine', :) = [];
asiT.relId = reshape(1:height(asiT), [], 1);

%% Build connectivity matrix
clear cur*;

pairT = table;
[pairT.ids, ~, curIds] = unique([ ...
    asiT.preAggloId, asiT.postAggloId], 'rows');

pairT.preAggloId = pairT.ids(:, 1);
pairT.postAggloId = pairT.ids(:, 2);
pairT(:, 'ids') = [];

pairT.axonClass = conn.axonMeta.axonClass(pairT.preAggloId);
pairT.targetClass = conn.denMeta.targetClass(pairT.postAggloId);

pairT.asiIds = accumarray( ...
    curIds, asiT.relId, [], @(ids) {ids});
pairT.meanLog10AsiArea = cellfun( ...
    @(ids) mean(log10(asiT.area(ids))), pairT.asiIds);

curCv = @(areas) std(areas) / mean(areas);
curMask = cellfun(@numel, pairT.asiIds) > 1;

pairT.cvAsiAreas(:) = nan;
pairT.cvAsiAreas(curMask) = cellfun( ...
    @(ids) curCv(asiT.area(ids)), pairT.asiIds(curMask));

%% Prepare for analysis
clear cur*;

curAxonClasses = setdiff( ...
    conn.axonMeta.axonClass, {'Inhibitory', 'Other'});
curTargetClasses = setdiff( ...
    conn.denMeta.targetClass, { ...
    'AxonInitialSegment', 'Ignore', 'Somata'});

curClassPairs = unique([ ...
    double(pairT.axonClass), ...
    double(pairT.targetClass)], 'rows');
curClassPairs = curClassPairs( ...
    ismember(curClassPairs(:, 1), double(curAxonClasses)) ...
  & ismember(curClassPairs(:, 2), double(curTargetClasses)), :);

plotConfigs = struct;
plotConfigs(1).classes = [ ...
    double(curAxonClasses(:)), ...
    zeros(size(curAxonClasses(:)))];

plotConfigs(2).classes = [ ...
    zeros(size(curTargetClasses(:))), ...
    double(curTargetClasses(:))];

plotConfigs(3).classes = curClassPairs;

%% Look synapse size vs. axon / dendrite classes and combinations thereof
clear cur*;

curBinEdges = linspace(-2, 0.5, 21);
for curPlotConfig = plotConfigs
    curClasses = curPlotConfig.classes;
    
    curVals = cell(size(curClasses, 1), 1);
    for curClassIdx = 1:size(curClasses, 1)
        curPreId = curClasses(curClassIdx, 1);
        curPostId = curClasses(curClassIdx, 2);
        
        curMask = ...
            (~curPreId | double(pairT.axonClass) == curPreId) ...
          & (~curPostId | double(pairT.targetClass) == curPostId);
        
        curVals{curClassIdx} = pairT.meanLog10AsiArea(curMask);
    end
    
    curAxonLeg = [{''}; categories(conn.axonMeta.axonClass)];
    curAxonLeg = curAxonLeg(1 + curClasses(:, 1));
    
    curDendLeg = [{''}; categories(conn.denMeta.targetClass)];
    curDendLeg = curDendLeg(1 + curClasses(:, 2));
    
    curLegends = arrayfun( ...
        @(ax, dend) strjoin(cat(2, ax, dend), ' → '), ...
        curAxonLeg, curDendLeg, 'UniformOutput', false);
    
    % Histogram
    curHistFig = figure();
    curHistFig.Color = 'white';
    
    curHistAx = axes(curHistFig); %#ok
    hold(curHistAx, 'on');
    
    cellfun(...
        @(data) histogram( ...
            curHistAx, data, ...
            'BinEdges', curBinEdges, ...
            'Normalization', 'probability', ...
            'DisplayStyle', 'stairs', 'LineWidth', 2), ...
        curVals);
    
    curHistLeg = legend(curHistAx, curLegends, 'Location', 'NorthWest');
    curHistLeg.Box = 'off';
    
    curHistAx.XLim = curBinEdges([1, end]);
    curHistAx.TickDir = 'out';
    
    % Boxplot
    curBoxFig = figure();
    curBoxFig.Color = 'white';
    
    curBoxAx = axes(curBoxFig); %#ok
    hold(curBoxAx, 'on');
    
    curGroupId = repelem(1:numel(curVals), cellfun(@numel, curVals));
    boxplot(cell2mat(curVals), curGroupId(:));
end

% Ratio of means
curT = table;
curT.title = curLegends;
curT.median = 10 .^ cellfun(@median, curVals);
curT.otherCcRatio = curT.median ./ ...
    curT.median(all(curClasses == [1, 4], 2));
curT.otherRatio = curT.median ./ arrayfun(@(a) ...
    curT.median(all(curClasses == [a, 4], 2)), curClasses(:, 1));
curT %#ok
