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

%% Look at remaining synapses from axons with small / large SASD pairs
clear cur*;

curClassName = 'Corticocortical';
curBinEdges = linspace(-1.5, 0.5, 21);

curRois = struct;
curRois(1).title = ...
    'Remaining synapses from axons with large area, low CV pair';
curRois(1).cvRange = [0, 0.4];
curRois(1).log10AsiRange = [-0.2, +0.2];

curRois(2).title = ...
    'Remaining synapses from axons with small area, low CV pair';
curRois(2).cvRange = [0.1, 0.5];
curRois(2).log10AsiRange = [-0.8, -0.5];

% Find bi-spiny connections within CV-log10(ASI) regions of interest
curPairIds = arrayfun( ...
    @(roi) find( ...
        pairT.cvAsiAreas >= roi.cvRange(1) ...
      & pairT.cvAsiAreas <= roi.cvRange(end) ...
      & pairT.meanLog10AsiArea >= roi.log10AsiRange(1) ...
      & pairT.meanLog10AsiArea <= roi.log10AsiRange(end) ...
      & cellfun(@numel, pairT.asiIds) == 2 ...
      & pairT.axonClass == curClassName), ...
	curRois, 'UniformOutput', false);
[curRois.pairIds] = deal(curPairIds{:});

curAxonIds = arrayfun( ...
    @(roi) unique(pairT.preAggloId(roi.pairIds)), ...
    curRois, 'UniformOutput', false);
[curRois.axonIds] = deal(curAxonIds{:});

curOtherAsiIds = arrayfun( ...
    @(roi) unique(cell2mat(pairT.asiIds( ...
        setdiff(find( ...
        ismember(pairT.preAggloId, roi.axonIds) ...
      & cellfun(@numel, pairT.asiIds) ~= 2), roi.pairIds)))), ...
    curRois, 'UniformOutput', false);

% Add fake ROI with all synapses
curRois(end + 1).title = 'All synapses';
curOtherAsiIds{end + 1} = unique(cell2mat(pairT.asiIds));

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [560, 545];
curAx = axes(curFig);
hold(curAx, 'on');

for curIdx = 1:numel(curRois)
    histogram(curAx, ...
        log10(asiT.area(curOtherAsiIds{curIdx})), ...
        'BinEdges', curBinEdges, 'Normalization', 'probability', ...
        'DisplayStyle', 'stairs', 'LineWidth', 2);
end

curLeg = legend(curAx, {curRois.title}, 'Location', 'SouthOutside');
curLeg.Box = 'off';

curAx.TickDir = 'out';
curAx.XLim = curBinEdges([1, end]);

xlim(curAx, curBinEdges([1, end]));
xlabel(curAx, 'Average log_{10}(ASI area [µm²])');
ylabel(curAx, 'Probability');

title(curAx, ...
    {info.filename; info.git_repos{1}.hash; sprintf( ...
    'Primary spine synapses from %s axons', curClassName)}, ...
    'FontSize', 10, 'FontWeight', 'normal');

% Report some numbers
for curIdx = 1:(numel(curRois) - 1)
    curRoi = curRois(curIdx);
    
   [~, curPVal] = kstest2( ...
       log10(asiT.area(curOtherAsiIds{curIdx})), ...
       log10(asiT.area(curOtherAsiIds{end})));
   
   fprintf('# %s\n', curRoi.title);
   fprintf('* CV range: %g to %g\n', curRoi.cvRange);
   fprintf('* Average log10(ASI) range: %g to %g\n', curRoi.cvRange);
   fprintf('* No. SASD pairs meeting criteria: %d\n', numel(curRoi.pairIds));
   fprintf('* No. axons with at least one such SASD pair: %d\n', numel(curRoi.axonIds));
   fprintf('* No. remaining synapses from above axons: %d\n', numel(curOtherAsiIds{curIdx}));
   fprintf('* Prob. of remaining synapses being different from overall: %g\n', curPVal);
   fprintf('\n');
end

%% Prototyping
clear cur*;

curCvRange = [0, 0.4];
curClassName = 'Corticocortical';

curPairT = pairT;
curPairT = curPairT(curPairT.axonClass == curClassName, :);
curPairT = sortrows(curPairT, 'meanLog10AsiArea', 'ascend');

curPairIds = find(...
    curPairT.cvAsiAreas >= curCvRange(1) ...
  & curPairT.cvAsiAreas <= curCvRange(end) ...
  & cellfun(@numel, curPairT.asiIds) == 2);

% Collected pairs
curCollPairIds = zeros(0, 1);
curMeanLog10Asi = zeros(0, 1);
curMeanLog10AsiPVal = zeros(0, 1);
curCvPVal = zeros(0, 1);

curFst = @(vals) vals(1);
curCv = @(areas) std(areas, 0, 2) ./ mean(areas, 2);

tic;
for curIdx = 1:numel(curPairIds)
    Util.progressBar(curIdx, numel(curPairIds));
    
    curPairId = setdiff(curPairIds, curCollPairIds);
    if isempty(curPairId); break; end
    curPairId = curPairId(1);
    
    curCollPairIds = cat(1, curCollPairIds, curPairId);
    
    % Complete
    curCollAxonIds = unique( ...
        curPairT.preAggloId(curCollPairIds));
    curCollPairIds = curPairIds(ismember( ...
        curPairT.preAggloId(curPairIds), curCollAxonIds));
    
    %{
    curCollAsiIds = unique(cell2mat(curPairT.asiIds( ...
        ismember(curPairT.preAggloId, curCollAxonIds))));
    %}
    
    curCollAsiIds = unique(cell2mat(curPairT.asiIds(curCollPairIds)));
    
    rng(0);
    curRandPairs = randi(numel(curCollAsiIds), 1000);
    curRandPairs = reshape(curRandPairs, [], 2);
    
    curRandCvs = curCv(asiT.area(curRandPairs));
    curRandMeanLog10AsiArea = mean(log10(asiT.area(curRandPairs)), 2);
    
    curRandMask = ...
        curRandCvs >= curCvRange(1) ...
      & curRandCvs <= curCvRange(end);
  
    curA = curRandMeanLog10AsiArea(curRandMask);
    curB = curPairT.meanLog10AsiArea(curCollPairIds);
   [~, curMeanLog10AsiPVal(curIdx)] = kstest2(curA, curB);
   
    curA = curRandCvs;
    curB = curPairT.cvAsiAreas(curCollPairIds);
   [~, curCvPVal(curIdx)] = kstest2(curA, curB);
    
    curMeanLog10Asi(curIdx) = curPairT.meanLog10AsiArea(curPairId);
end

% Plotting
curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [560, 540];

curAx = subplot(3, 1, 1);
plot(curAx, log10(curCvPVal));
ylabel(curAx, {'log(p-value)'; 'for CV distribution'});

curAx = subplot(3, 1, 2);
plot(curAx, log10(curMeanLog10AsiPVal));
ylabel(curAx, {'log(p-value)'; 'for low-CV ASI areas'});

curAx = subplot(3, 1, 3);
plot(curAx, curMeanLog10Asi);
ylabel(curAx, {'Average log_{10}(ASI)'; 'of SASD pair'});
xlabel(curAx, 'SASD size threshold / axons');

set(curFig.Children, 'TickDir', 'out');

annotation( ...
    curFig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
