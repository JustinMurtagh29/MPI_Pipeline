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

%% Build plot configurations
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

%% Look at synapse sizes
clear cur*;

for curPlotConfig = plotConfigs
    curClasses = curPlotConfig.classes;
    
    % Collect axon-spine interface areas
    curAsiIds = cell(size(curClasses, 1), 1);
    for curClassIdx = 1:size(curClasses, 1)
        curPreId = curClasses(curClassIdx, 1);
        curPostId = curClasses(curClassIdx, 2);
        
        curMask = ...
            (~curPreId | double(pairT.axonClass) == curPreId) ...
          & (~curPostId | double(pairT.targetClass) == curPostId);
        curAsiIds{curClassIdx} = unique(cell2mat(pairT.asiIds(curMask)));
    end
    
    % Build legends
    curAxonLeg = [{''}; categories(conn.axonMeta.axonClass)];
    curAxonLeg = curAxonLeg(1 + curClasses(:, 1));
    
    curDendLeg = [{''}; categories(conn.denMeta.targetClass)];
    curDendLeg = curDendLeg(1 + curClasses(:, 2));
    
    curLegends = cellfun( ...
        @(ax, dend, n) sprintf('%s → %s (n = %d)', ax, dend, n), ...
        curAxonLeg, curDendLeg, num2cell(cellfun(@numel, curAsiIds)), ...
        'UniformOutput', false);
    
    % Boxplot
    curBoxFig = figure();
    curBoxFig.Color = 'white';
    
    curBoxAx = axes(curBoxFig); %#ok
    hold(curBoxAx, 'on');
    
    curVals = log10(asiT.area(cell2mat(curAsiIds)));
    curGroupIds = repelem(1:numel(curAsiIds), cellfun(@numel, curAsiIds));
    boxplot(curVals(:), curGroupIds(:));
    
    curBoxAx.Box = 'off';
    curBoxAx.TickDir = 'out';
    
    ylabel(curBoxAx, 'log_{10}(axon-spine interface area [µm²])');
    xticklabels(curBoxAx, curLegends);
    curBoxAx.XTickLabelRotation = 20;
    
    title(curBoxAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontSize', 10, 'FontWeight', 'normal');
    
    % Quantitative evaluation
    curPVals = nan(size(curClasses, 1));
    curMedRatios = nan(size(curClasses, 1));
    for curIdxA = 1:size(curClasses, 1)
        for curIdxB = curIdxA:size(curClasses, 1)
            
            curAreasA = asiT.area(curAsiIds{curIdxA});
            curAreasB = asiT.area(curAsiIds{curIdxB});
            
           [~, curPVal] = kstest2(curAreasA, curAreasB);
            curPVals(curIdxA, curIdxB) = curPVal;
            curPVals(curIdxB, curIdxA) = curPVal;
           
            curMedAreaA = median(curAreasA);
            curMedAreaB = median(curAreasB);
            
            curMedRatios(curIdxA, curIdxB) = curMedAreaA / curMedAreaB;
            curMedRatios(curIdxB, curIdxA) = curMedAreaB / curMedAreaA;
        end
    end
    
    curPVals = array2table(curPVals);
    curPVals.Properties.RowNames = curLegends;
    curPVals %#ok
    
    curMedRatios = array2table(curMedRatios);
    curMedRatios.Properties.RowNames = curLegends;
    curMedRatios %#ok
end

%% Look at synapse size similarity
clear cur*;

% Cosmetics
curBinEdges = linspace(0, 1.5, 16);

for curPlotConfig = plotConfigs
    curClasses = curPlotConfig.classes;
    
    % Collect axon-spine interface areas
    curAsiIds = cell(size(curClasses, 1), 1);
    curSasdIds = cell(size(curClasses, 1), 1);
    curCtrlIds = cell(size(curClasses, 1), 1);
    
    for curClassIdx = 1:size(curClasses, 1)
        curPreId = curClasses(curClassIdx, 1);
        curPostId = curClasses(curClassIdx, 2);
        
        curMask = ...
            (~curPreId | double(pairT.axonClass) == curPreId) ...
          & (~curPostId | double(pairT.targetClass) == curPostId);
        curAsiIds{curClassIdx} = unique(cell2mat(pairT.asiIds(curMask)));
        
        curMask = curMask & cellfun(@numel, pairT.asiIds) == 2;
        curSasd = cell2mat(pairT.asiIds(curMask));
        curSasd = transpose(reshape(curSasd, 2, []));
        curSasdIds{curClassIdx} = curSasd;
    
        % Build random pairs
        rng(0);
        curRands = curAsiIds{curClassIdx};
        curRands = curRands(randperm(numel(curRands)));
        curRands = curRands(1:(2 * floor(numel(curRands) / 2)));
        curRands = reshape(curRands, [], 2);
        curCtrlIds{curClassIdx} = curRands;
    end
    
    % Remove empty classes
    curMask = cellfun(@isempty, curSasdIds);
    curClasses(curMask, :) = [];
    curAsiIds(curMask, :) = [];
    curSasdIds(curMask, :) = [];
    curCtrlIds(curMask, :) = [];
    
    % Build legends
    curAxonLeg = [{''}; categories(conn.axonMeta.axonClass)];
    curAxonLeg = curAxonLeg(1 + curClasses(:, 1));
    
    curDendLeg = [{''}; categories(conn.denMeta.targetClass)];
    curDendLeg = curDendLeg(1 + curClasses(:, 2));
    
    curLegends = cellfun( ...
        @(ax, dend, n) sprintf('%s → %s (n = %d)', ax, dend, n), ...
        curAxonLeg, curDendLeg, num2cell(cellfun(@numel, curAsiIds)), ...
        'UniformOutput', false);
    
    curCvs = @(vals) ...
        std(vals, 0, 2) ./ mean(vals, 2);
    curCvs = cellfun( ...
        @(ids) curCvs(asiT.area(ids)), ...
        [curSasdIds, curCtrlIds], 'UniformOutput', false);
    
   [~, curPVals] = cellfun( ...
       @(a, b) kstest2(a, b, 'tail', 'larger'), ...
       curCvs(:, 1), curCvs(:, 2), 'ErrorHandler', ...
       @(~, ~, ~) deal(false, nan));
   [curLearnedFrac, curUnlearnedFrac, curCvThresh] = cellfun( ...
        @(a, b) connectEM.Consistency.calculateLearnedFraction( ...
            asiT, struct('synIdPairs', a), struct('synIdPairs', b), ...
            'method', 'maxdifference'), ...
        curSasdIds, curCtrlIds);
    
    curRes = { ...
        'pVal', 'numPairs', 'cvTresh', 'learnedFrac', 'unlearnedFrac'};
    curRes = array2table( ...
       [curPVals, cellfun(@numel, curSasdIds) / 2, ...
        curCvThresh, curLearnedFrac, curUnlearnedFrac], ...
        'VariableNames', curRes, 'RowNames', curLegends) %#ok
    
    for curRow = 1:size(curCvs, 1)
        curFig = figure;
        curFig.Color = 'white';
        curFig.Position(3:4) = [540, 570];
        
        curAx = axes(curFig); %#ok
        hold(curAx, 'on');
        
        cellfun( ...
            @(vals) histogram( ...
                curAx, vals, ...
                'BinEdges', curBinEdges, ...
                'Normalization', 'probability', ...
                'DisplayStyle', 'stairs', ...
                'LineWidth', 2), ...
            curCvs(curRow, :));
        
        curAx.TickDir = 'out';
        curAx.XLim = curBinEdges([1, end]);
        axis(curAx, 'square');
        xlabel(curAx, 'Coefficient of variation (CV)');
        ylabel(curAx, 'Probability');
        
        curLeg = { ...
            sprintf( ...
                'Same-axon same-dendrite pairs (n = %d)', ...
                size(curSasdIds{curRow}, 1)), ...
            sprintf( ...
                'Random pairs (n = %d)', ...
                size(curCtrlIds{curRow}, 1))};
        curLeg = legend(curAx, curLeg, 'Location', 'NorthEast');
        curLeg.Box = 'off';
        
        curTitle = sprintf( ...
           ['p = %.3g, CV = %.2g, ', ...
            'learned = %.1f %%, unlearned = %.1f %%'], ...
            curPVals(curRow), curCvThresh(curRow), ...
            100 * [curLearnedFrac(curRow), curUnlearnedFrac(curRow)]);
        curTitle = { ...
            info.filename; info.git_repos{1}.hash; ...
            curLegends{curRow}; curTitle}; %#ok
        title(curAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
    end
    
    curTitles = cellfun(@(legend, asiIds) regexprep( ...
        legend, 'n = \d+', sprintf('n = %d pairs', numel(asiIds) / 2)), ...
        curLegends, curSasdIds, 'UniformOutput', false);
    curShuffle = @(vals) vals(reshape(randperm(numel(vals)), [], 2));
    cellfun(@(title, asiIds) ...
        connectEM.Consistency.plotSizeSimilarityHeatmap( ...
            asiT, ...
            struct('synIdPairs', asiIds), ...
            struct('synIdPairs', curShuffle(asiIds)), ...
            'title', {info.filename; info.git_repos{1}.hash; title}), ...
        curTitles, curSasdIds);
end

%% Actually analyse weight matrix
clear cur*;

%{
curAxonClass = 'Corticocortical';
curTargetClass = 'OtherDendrite';

curPairT = pairT( ...
    pairT.axonClass == curAxonClass ...
  & pairT.targetClass == curTargetClass, :);

[curAxonIds, ~, curPairT.preAggloId] = unique(curPairT.preAggloId);
[curDendIds, ~, curPairT.postAggloId] = unique(curPairT.postAggloId);

curPairT.medLog10AsiArea = cellfun(@(ids) ...
    median(log10(asiT.area(ids))), curPairT.asiIds);
%}

%% Generate fake data
clear cur*;

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

%%
curAreas = cellfun(@median, curConn);
curMu = mean(curAreas(:), 'omitnan');
curM = curAreas - curMu;
curM(isnan(curM)) = 0;

curNorm = 1 - isnan(curAreas);
curNorm = curNorm * curNorm';
curCorr = (curM * curM') ./ curNorm;
curCorr(isinf(curCorr)) = nan;

curLink = linkage(curCorr, 'average');
curSort = connectEM.Consistency.linkageToSorting(curLink);
figure; imagesc(curCorr(curSort, curSort)); axis square

%% Look at remaining synapses from axons with small / large SASD pairs
%{
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
    curCollAsiIds = unique(cell2mat( ...
        curPairT.asiIds(curCollPairIds)));
    
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
%}
