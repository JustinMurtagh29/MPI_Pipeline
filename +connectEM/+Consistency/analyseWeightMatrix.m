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
    connectEM.Consistency.loadConnectome(param, 'specificityClasses');

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
pairT.medianLog10AsiArea = cellfun( ...
    @(ids) median(log10(asiT.area(ids))), pairT.asiIds);

curCv = @(areas) std(areas) / mean(areas);
curMask = cellfun(@numel, pairT.asiIds) > 1;

pairT.cvAsiAreas(:) = nan;
pairT.cvAsiAreas(curMask) = cellfun( ...
    @(ids) curCv(asiT.area(ids)), pairT.asiIds(curMask));

%% Prepare for analysis
clear cur*;

curAxonClasses = setdiff( ...
    conn.axonMeta.axonClass, {'Ignore', 'Other'});
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
    
    curCvs = cell(size(curClasses, 1), 1);
    for curClassIdx = 1:size(curClasses, 1)
        curPreId = curClasses(curClassIdx, 1);
        curPostId = curClasses(curClassIdx, 2);
        
        curMask = ...
            (~curPreId | double(pairT.axonClass) == curPreId) ...
          & (~curPostId | double(pairT.targetClass) == curPostId);
        
        curCvs{curClassIdx} = pairT.medianLog10AsiArea(curMask);
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
        curCvs);
    
    curHistLeg = legend(curHistAx, curLegends, 'Location', 'NorthWest');
    curHistLeg.Box = 'off';
    
    curHistAx.XLim = curBinEdges([1, end]);
    curHistAx.TickDir = 'out';
    
    xlabel(curHistAx, 'Median ASI area of connection');
    ylabel(curHistAx, 'Probability');
    
    title(curHistAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontSize', 10, 'FontWeight', 'normal');
    
    % Boxplot
    curBoxFig = figure();
    curBoxFig.Color = 'white';
    
    curBoxAx = axes(curBoxFig); %#ok
    hold(curBoxAx, 'on');
    
    curGroupId = repelem(1:numel(curCvs), cellfun(@numel, curCvs));
    boxplot(cell2mat(curCvs), curGroupId(:));
    
    curBoxAx.Box = 'off';
    curBoxAx.TickDir = 'out';
    
    ylabel(curBoxAx, 'Median ASI area of connection');
    xticklabels(curBoxAx, curLegends);
    curBoxAx.XTickLabelRotation = 20;
    
    title(curBoxAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontSize', 10, 'FontWeight', 'normal');
end

% Ratio of medians
curT = table;
curT.title = curLegends;
curT.n = cellfun(@numel, curCvs);
curT.median = 10 .^ cellfun(@median, curCvs);
curT.otherCcRatio = curT.median ./ ...
    curT.median(all(curClasses == [1, 4], 2));
curT.ccRatio = curT.median ./ arrayfun(@(a) ...
    curT.median(all(curClasses == [1, a], 2)), curClasses(:, 2));
curT.otherRatio = curT.median ./ arrayfun(@(a) ...
    curT.median(all(curClasses == [a, 4], 2)), curClasses(:, 1));
curT %#ok

% For comparison, ratio between median TC vs. CC synapse size
curMedTc = median(asiT.area( ...
    conn.axonMeta.axonClass(asiT.preAggloId) == 'Thalamocortical'));
curMedCc = median(asiT.area( ...
    conn.axonMeta.axonClass(asiT.preAggloId) == 'Corticocortical'));
curMedTcCcRatio = curMedTc / curMedCc %#ok

%% Look at synapse size similarity
clear cur*;

curPValThresh = 0.10;
curBinEdges = linspace(0, 1.5, 16);

for curPlotConfig = plotConfigs
    curClasses = curPlotConfig.classes;
    
    curCvs = cell(size(curClasses, 1), 1);
    curMedAsis = cell(size(curClasses, 1), 1);
    curNullCvs = cell(size(curClasses, 1), 1);
    curNullMedAsis = cell(size(curClasses, 1), 1);
    
    for curClassIdx = 1:size(curClasses, 1)
        curPreId = curClasses(curClassIdx, 1);
        curPostId = curClasses(curClassIdx, 2);
        
        curMask = ...
            cellfun(@numel, pairT.asiIds) == 2 ...
          & (~curPreId | double(pairT.axonClass) == curPreId) ...
          & (~curPostId | double(pairT.targetClass) == curPostId);
      
        curCvs{curClassIdx} = pairT.cvAsiAreas(curMask);
        curMedAsis{curClassIdx} = pairT.medianLog10AsiArea(curMask);
        
        % Build random pairs
        curAsis = asiT.area(cell2mat(pairT.asiIds(curMask)));
        
        if isempty(curAsis)
            curNullCvs{curClassIdx} = nan(0, 1);
            curNullMedAsis{curClassIdx} = nan(0, 1);
            continue;
        end
        
        rng(0);
        curRandPairs = randi(numel(curAsis), [10000, 2]);
        curRandPairs = curAsis(curRandPairs);
        
        curNullCvs{curClassIdx} = ...
            std(curRandPairs, 0, 2) ...
         ./ mean(curRandPairs, 2);
        curNullMedAsis{curClassIdx} = ...
            median(log10(curRandPairs), 2);
    end
    
    curAxonLeg = [{''}; categories(conn.axonMeta.axonClass)];
    curAxonLeg = curAxonLeg(1 + curClasses(:, 1));
    
    curDendLeg = [{''}; categories(conn.denMeta.targetClass)];
    curDendLeg = curDendLeg(1 + curClasses(:, 2));
    
    curLegends = arrayfun( ...
        @(ax, dend) strjoin(cat(2, ax, dend), ' → '), ...
        curAxonLeg, curDendLeg, 'UniformOutput', false);
    
    % NOTE(amotta): Let's check if CVs are significantly smaller than
    % expected. If yes (using a certain p-value threshold), we determine
    % the CV value at which the difference between the cumulative density
    % functions is largest. This corresponds to the intersection of the
    % probability density functions.
    
    curT = table;
    curT.name = curLegends;
    curT.n = cellfun(@numel, curCvs);
    
    % NOTE(amotta): `kstest2` throws an error if any of the inputs is
    % empty. Let's handle this case by returning NaN as p-value.
   [~, curT.pVal] = cellfun( ...
       @(obs, null) kstest2(obs, null, 'tail', 'larger'), ...
       curCvs, curNullCvs, 'ErrorHandler', @(varargin) deal(false, nan));
   
    curSumTwo = @(v) [v(:, 1), cumsum(v(:, 2))];
    curCdfDiffs = cellfun(@(obs, null) ...
        curSumTwo(sortrows([ ...
           [obs, repmat(+1 / numel(obs), size(obs))]; ...
           [null, repmat(-1 / numel(null), size(null))]], 1)), ...
            curCvs, curNullCvs, 'UniformOutput', false);
    
    curT.cvThresh = nan(size(curCdfDiffs));
    curT.cvThresh(curT.pVal < curPValThresh) = cellfun( ...
        @(v) v(Util.nthOutput(2, @() max(v(:, 2))), 1), ...
        curCdfDiffs(curT.pVal < curPValThresh));
    
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
        curCvs);
    
    curHistLeg = legend(curHistAx, curLegends, 'Location', 'NorthEast');
    curHistLeg.Box = 'off';
    
    curHistAx.XLim = curBinEdges([1, end]);
    curHistAx.TickDir = 'out';
    
    xlabel(curHistAx, 'Coefficient of variation');
    ylabel(curHistAx, 'Probability');
    
    title(curHistAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontSize', 10, 'FontWeight', 'normal');
    
    % Boxplot
    curBoxFig = figure();
    curBoxFig.Color = 'white';
    
    curBoxAx = axes(curBoxFig); %#ok
    hold(curBoxAx, 'on');
    
    curGroupId = repelem(1:numel(curCvs), cellfun(@numel, curCvs));
    boxplot(cell2mat(curCvs), curGroupId(:));
    
    curBoxAx.Box = 'off';
    curBoxAx.TickDir = 'out';
    
    ylabel(curBoxAx, 'Coefficient of variation');
    xticklabels(curBoxAx, curLegends(cellfun(@numel, curCvs) > 0));
    curBoxAx.XTickLabelRotation = 20;
    
    title(curBoxAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontSize', 10, 'FontWeight', 'normal');
        
    % Bar plot
    curT.nullFracBelowThresh = cellfun( ...
        @(cvs, t) mean(cvs < t), curNullCvs, num2cell(curT.cvThresh));
    curT.obsFracBelowThresh = cellfun( ...
        @(cvs, t) mean(cvs < t), curCvs, num2cell(curT.cvThresh));
    curT.diffFracBelowThresh = ...
        curT.obsFracBelowThresh - curT.nullFracBelowThresh;
    
   [~, curT.areaP] = cellfun( ...
       @(aO, cO, aN, cN, t) kstest2(aO(cO < t), aN(cN < t)), ...
       curMedAsis, curCvs, curNullMedAsis, ...
       curNullCvs, num2cell(curT.cvThresh), ...
       'ErrorHandler', @(varargin) deal(false, nan));
    curT %#ok
    
    %{
    % Plot size distribution in low-CV area
    for curIdx = reshape(find(curT.areaP < 0.1), 1, [])
        curBins = linspace(-2, 0.5, 26);
        figure;
        hold on;
        histogram( ...
            curMedAsis{curIdx}( ...
                curCvs{curIdx} < curT.cvThresh(curIdx)), ...
            'Normalization', 'probability', 'BinEdges', curBins, ...
            'DisplayStyle', 'stairs', 'LineWidth', 2);
        histogram( ...
            curNullMedAsis{curIdx}( ...
                curNullCvs{curIdx} < curT.cvThresh(curIdx)), ...
            'Normalization', 'probability', 'BinEdges', curBins, ...
            'DisplayStyle', 'stairs', 'LineWidth', 2);
    end
    %}
    
    curBarFig = figure();
    curBarFig.Color = 'white';
    
    curBarAx = axes(curBarFig); %#ok
    hold(curBarAx, 'on');
    
    curColors = get(groot, 'defaultAxesColorOrder');
    bar(curBarAx, curT.nullFracBelowThresh, 'FaceColor', 'none', ...
        'EdgeColor', 'black', 'LineWidth', 2);
    bar(curBarAx, curT.obsFracBelowThresh, 'FaceColor', 'none', ...
        'EdgeColor', curColors(1, :), 'LineWidth', 2);
    
    curBarAx.Box = 'off';
    curBarAx.TickDir = 'out';
    
    ylim(curBarAx, [0, 1]);
    xlim(curBarAx, [0, numel(curCvs) + 1]);
    xticks(curBarAx, 1:numel(curCvs));
    xticklabels(curBarAx, curLegends);
    curBarAx.XTickLabelRotation = 20;
    
    ylabel(curBarAx, { ...
        'Fraction of bi-synaptic connections'; ...
        'with unexpetedly small size variability'});
    
    curLeg = {'Null model', 'Observed'};
    curLeg = legend(curBarAx, curLeg, 'Location', 'NorthEast');
    curLeg.Box = 'off';
    
    title(curBarAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontSize', 10, 'FontWeight', 'normal');
end
