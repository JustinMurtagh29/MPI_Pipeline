% This is a minimally modified version of
% +connectEM/+Connectome/plotSynapseSizeConsistency.m
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
connFile = fullfile(rootDir, 'connectome', 'Connectome_20191203T021242-results_20191203T021242-results-auto-spines-v2_SynapseAgglomerates--20191203T021242-results--20191203T021242-results-auto-spines-v2--v1.mat');
asiRunId = '20191213T143516';

% NOTE(amotta): Set this path to generate a debug NML file
debugNmlDir = '/home/amotta/Desktop/';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = load(connFile);

synFile = conn.info.param.synFile;
syn = load(synFile);

curAsiFile = sprintf('__%s_asiT.mat', asiRunId);
curAsiFile = strrep(connFile, '.mat', curAsiFile);
asiT = Util.load(curAsiFile, 'asiT');

% NOTE(amotta): From +connectEM/+Connectome/plotSynapseSizeConsistency.m
asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

curTypeFile = strrep(synFile, '.mat', '__types_v1.mat');
syn.synapses.type = Util.load(curTypeFile, 'types');

%% Build synapse table
clear cur*;

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

%% Build axon classes
clear cur*;

curMeta = conn.axonMeta;
% NOTE(amotta): Remove "old" synapse counts
curVars = curMeta.Properties.VariableNames;
curVars = setdiff(curVars, {'synCount', 'spineSynCount'}, 'stable');
curMeta = curMeta(:, curVars);

curMeta.synCount = accumarray( ...
    synT.preAggloId, 1, [height(curMeta), 1]);
curMeta.priSpineSynCount = accumarray( ...
    synT.preAggloId, synT.type == 'PrimarySpine', [height(curMeta), 1]);
curMeta.priSpineSynFrac = curMeta.priSpineSynCount ./ curMeta.synCount;

curMeta.class(:) = categorical({'Unknown'});
curSpeckle = curMeta.synCount < 10;
curMeta.class(curSpeckle) = 'Speckle';
curMergerMask = curMeta.synCount > 25;
curMeta.class(curMergerMask) = 'Merger';
curOkayMask = not(curSpeckle | curMergerMask);
curMeta.class(curOkayMask & curMeta.priSpineSynFrac > 0.5) = 'Excitatory';
curMeta.class(curOkayMask & curMeta.priSpineSynFrac < 0.25) = 'Inhibitory';

conn.axonMeta = curMeta;

%% Complete ASI table
clear cur*;

[~, curTypeRow] = ismember(asiT.id, synT.id);
asiT.type = synT.type(curTypeRow);
asiT.axonClass = conn.axonMeta.class(asiT.preAggloId);

%% Build configurations
clear cur*;

plotConfigs = struct('synIds', {}, 'title', {}, 'tag', {});

axonClasses = {'Exc', {'Excitatory'}};
targetClasses = {'All', 'OtherDendrite'};

for curAxonIdx = 1:size(axonClasses, 1)
    curAxonClass = axonClasses(curAxonIdx, :);
    
    for curTargetIdx = 1:size(targetClasses, 1)
        curTargetClass = targetClasses(curTargetIdx, :);
        
        curSynIds = find( ...
            asiT.type == 'PrimarySpine' ...
          & ismember(asiT.axonClass, curAxonClass{2}) ...
          & ismember(asiT.targetClass, curTargetClass{2}));
        if isempty(curSynIds); continue; end
        
        curTitle = sprintf( ...
            '%s → %s primary spine synapses', ...
            curAxonClass{1}, curTargetClass{1});
        curTag = sprintf('%s %s pri sp', ...
            curAxonClass{1}, curTargetClass{1});
        
        curPlotConfig = plotConfigs([]);
        curPlotConfig(1).synIds = curSynIds;
        curPlotConfig(1).title = curTitle;
        curPlotConfig(1).tag = curTag;
        
        plotConfigs(end + 1) = curPlotConfig; %#ok
    end
end

plotConfigs = reshape( ...
    plotConfigs, ...
    size(targetClasses, 1), ...
    size(axonClasses, 1));

%% Set p-value thresholds for evaluation
clear cur*;

% NOTE(amotta): p-value threshold at which large and small low-CV regions
% are not yet merged. (I've only looked at range up to p = 10 %.)
curPvalThreshs = [0.5, 1, 2, 3, 4, 5] / 100;
plotConfigs(1, 1).twoDimPValThreshs = 0.0102;

for curIdx = 1:numel(plotConfigs)
    curThreshs = plotConfigs(curIdx).twoDimPValThreshs;
    
    % Only consider p-value thresholds before LTP and LTD regions merger
    curThreshs = curPvalThreshs(curPvalThreshs < curThreshs);
    
    % Only consider lowest and highest p-value threshold
    curThreshs = curThreshs(unique([1, numel(curThreshs)]));
    plotConfigs(curIdx).twoDimPValThreshs = curThreshs;
end

%% Export random examples to webKnossos
clear cur*;
rng(0);

curPairCount = 50;
curMaxNodeCount = 50000;

if ~isempty(debugNmlDir)
    if ~exist('axons', 'var')
        axons = conn.info.param.axonFile;
        axons = Util.load(axons, 'axons');
    end
    
    if ~exist('dendrites', 'var')
        dendrites = conn.info.param.dendriteFile;
        dendrites = Util.load(dendrites, 'dendrites');
    end
    
    curPlotConfig = plotConfigs(1);
    curPairs = connectEM.Consistency.buildPairConfigs(asiT, curPlotConfig);
    
    curPairT = table;
    curPairT.synIds = curPairs(1).synIdPairs;
    curPairT.axonId = asiT.preAggloId(curPairT.synIds(:, 1));
    curPairT.dendId = asiT.postAggloId(curPairT.synIds(:, 1));
    
    curPairCount = min(curPairCount, height(curPairT));
    curRandIds = randperm(height(curPairT), curPairCount);
    curNumDigits = ceil(log10(1 + curPairCount));
    
    for curIdx = 1:curPairCount
        curId = curRandIds(curIdx);
        curAxonId = curPairT.axonId(curId);
        curDendId = curPairT.dendId(curId);
        
        curAsiIds = curPairT.synIds(curId, :);
        curAsiPos = asiT.pos(curAsiIds, :);
        
        curAxon = axons(curAxonId);
        curDend = dendrites(curDendId);
        
        if size(curAxon.nodes, 1) > curMaxNodeCount ...
                || size(curDend.nodes, 1) > curMaxNodeCount
            warning('Pair %d skipped due to agglomerate size', curIdx);
            continue;
        end
        
        curSkel = skeleton();
        curSkel = Skeleton.setParams4Pipeline(curSkel, param);
        curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
        
        curSkel = curSkel.addTree( ...
            sprintf('Axon %d', curAxonId), ...
            curAxon.nodes(:, 1:3), curAxon.edges);
        curSkel = curSkel.addTree( ...
            sprintf('Dendrite %d', curDendId), ...
            curDend.nodes(:, 1:3), curDend.edges);
        
        for curAsiIdx = 1:numel(curAsiIds)
            curSkel = curSkel.addTree( ...
                sprintf('ASI %d', curAsiIds(curAsiIdx)), ...
                curAsiPos(curAsiIdx, :));
        end
        
        curNmlFile = fullfile(debugNmlDir, sprintf( ...
            '%0*d_sasd-pair-%d.nml', curNumDigits, curIdx, curId));
        curSkel.write(curNmlFile);
    end
end

%% Synapse sizes
clear cur*;

connectEM.Consistency.plotSizeHistogram( ...
    info, asiT, plotConfigs, 'scale', 'log10');

%% Report number of occurences for degree of coupling
clear cur*;

curPlotConfig = plotConfigs(1, 1);
curSynT = asiT(curPlotConfig.synIds, :);
curDoc = curSynT(:, {'preAggloId', 'postAggloId'});
[~, ~, curDoc] = unique(curDoc, 'rows');
curDoc = accumarray(curDoc, 1);

docT = table;
[docT.degreeOfCoupling, ~, curDocCount] = unique(curDoc);
docT.occurences = accumarray(curDocCount, 1);
docT(docT.degreeOfCoupling == 1, :) = [];

fprintf('Degree of coupling histogram\n\n');
disp(docT)

%% Plot histogram of degree of coupling
connectEM.Consistency.plotCouplingHistogram( ...
    info, asiT, plotConfigs(1, 1), 'normalization', 'count');

%% Synapse area variability
clear cur*;
pValues = nan(4, numel(plotConfigs));

for curIdx = 1:numel(plotConfigs)
    curPlotConfig = plotConfigs(curIdx);
    
    curPairConfigs = ...
        connectEM.Consistency.buildPairConfigs(asiT, curPlotConfig);
    
    curFig = ...
        connectEM.Consistency.plotVariabilityHistogram( ...
            info, asiT, curPlotConfig, curPairConfigs(:));
    curFig.Position(3:4) = [370, 540];
   [curLearnedFrac, curUnlearnedFrac, curCvThresh] = ...
        connectEM.Consistency.calculateLearnedFraction( ...
            asiT, curPairConfigs(1), curPairConfigs(end), ...
            'method', 'maxdifference');
    curOpenFrac = 1 - curLearnedFrac - curUnlearnedFrac;

    % Quantitative reporting
    fprintf('**%s**\n', curPlotConfig.title);
    fprintf('CV between pairs (mean ± std)\n');
    
    for curPairConfig = curPairConfigs
        curCvs = asiT.area(curPairConfig.synIdPairs);
        curCvs = std(curCvs, 0, 2) ./ mean(curCvs, 2);

        fprintf( ...
            '* %s: %f ± %f\n', ...
            curPairConfig.title, ...
            mean(curCvs), std(curCvs));
    end
    
    fprintf('\n');
    
    fprintf('Areas under and between curves\n');
    fprintf('* Learned fraction: %.1f %%\n', 100 * curLearnedFrac);
    fprintf('* Possibly learned fraction: %.1f %%\n', 100 * curOpenFrac);
    fprintf('* Unlearned fraction: %.1f %%\n', 100 * curUnlearnedFrac);
    fprintf('* CV threshold: %.2f\n', curCvThresh);
    fprintf('\n');

    fprintf('Significance tests\n');
    fprintf('p-values for unexpected synapse size similarity\n');
    
    for curPairIdx = 1:4
        curPairConfig = curPairConfigs(curPairIdx);
        
        curPValue = ...
            connectEM.Consistency.testVariability( ...
                asiT, curPairConfig, curPairConfigs(end));
        pValues(curPairIdx, curIdx) = curPValue;
        
        fprintf('* %s: %g\n', curPairConfig.title, curPValue);
    end

    fprintf('\n');
    fprintf('\n');
end

pValues = reshape(pValues, [4, size(plotConfigs)]);

%% Correlation of synapse size correlation with synapse size
clear cur*;

% Density difference map
curAxisX = 'cv';
curScaleY = 'log10';
curImSize = [256, 256];
curMethod = 'kde2d';

% Show region in which the same-axon same-dendrite primary spine synapse
% pairs with at least one synapse in the smallest 10th percentile are
% located.
curMinPrctiles = 10;

switch lower(curAxisX)
    case 'cv'
        curLimX = [0, 1.5];
        curTicksX = linspace(0, 1.5, 4);
        curFuncX = @(a) std(a, 0, 2) ./ mean(a, 2);
        curLabelX = 'Coefficient of variation';
    case 'reldiff'
        curLimX = [0, 2];
        curTicksX = linspace(0, 2, 5);
        curFuncX = @(a) abs(diff(a, 1, 2)) ./ mean(a, 2);
        curLabelX = 'Relative difference';
    otherwise
        error('Invalid X axis "%s"', curAxisX);
end

switch lower(curScaleY)
    case 'linear'
        curLimY = [0, 1];
        curTicksY = linspace(0, 1, 5);
        curFuncY = @(areas) mean(areas, 2);
        curLabelY = 'Average ASI area [µm²]';
    case 'log10'
        curLimY = [-1.5, 0.5];
        curTicksY = linspace(-1.5, 0.5, 5);
        curFuncY = @(areas) log10(mean(areas, 2));
        curLabelY = 'log10(Average ASI area [µm²])';
    otherwise
        error('Invalid Y scale "%s"', curScaleY);
end

%{
% NOTE(amotta): The following color scheme was use for the scatter plot
% shown in presentation at the Connectomics Conference 2019 in Berlin.
%
% Region 1: Low CV and small. Red.
% Region 2: Low CV and large. Green.
curReds = autumn(2 + 2);
curGreens = summer(2 + 2);

curRegionColors = cell(1, 2);
curRegionColors{1} = curReds([1, end - 1], :);
curRegionColors{2} = curGreens([1, end - 1], :);
%}

curRegionColors = [];
curRegionContourProps = {'LineColor', 'black'};

% NOTE(amotta): The `curMinMap` and `curMaxMap` matrices have the same size
% as the heatmaps of the CV × log10(avg. ASI) space. They contain the ASI
% areas of the smaller and larger synapses, respectively.
curLog10Avg = linspace(curLimY(1), curLimY(2), curImSize(1));
if strcmpi(curScaleY, 'linear'); curLog10Avg = log10(curLog10Avg); end
curCv = linspace(curLimX(1), curLimX(2), curImSize(2));
if strcmpi(curAxisX, 'reldiff'); curCv = curCv / sqrt(2); end
curCv(curCv >= sqrt(2)) = nan;

curMinMap = (10 .^ curLog10Avg(:)) .* (1 - curCv / sqrt(2));
curMaxMap = (10 .^ curLog10Avg(:)) .* (1 + curCv / sqrt(2));

for curConfig = reshape(plotConfigs(1, :), 1, [])
    curSaSdConfig = connectEM.Consistency.buildPairConfigs(asiT, curConfig);
    curSaSdConfig = curSaSdConfig(1);
            
    % Connection types
    curSaSdT = table;
    curSaSdT.synIds = curSaSdConfig.synIdPairs;
    curSaSdT.axonClass = asiT.axonClass(curSaSdT.synIds(:, 1));
    curSaSdT.targetClass = asiT.targetClass(curSaSdT.synIds(:, 1));

    curSaSdT.areas = asiT.area(curSaSdT.synIds);
    curSaSdT.x = curFuncX(curSaSdT.areas);
    curSaSdT.y = curFuncY(curSaSdT.areas);
    
    curSaSdT = curSaSdT( ...
        curLimX(1) <= curSaSdT.x & curSaSdT.x <= curLimX(2) ...
      & curLimY(1) <= curSaSdT.y & curSaSdT.y <= curLimY(2), :);
    
    curSaSdT.mapIds = [ ...
        discretize(curSaSdT.x, linspace( ...
            curLimX(1), curLimX(2), curImSize(2) + 1)), ...
        discretize(curSaSdT.y, linspace( ...
            curLimY(1), curLimY(2), curImSize(1) + 1))];
    curSaSdT.mapIdx = sub2ind( ...
        curImSize, curSaSdT.mapIds(:, 2), curSaSdT.mapIds(:, 1));
    
    % Restrict analysis to convex hull around point cloud
    curRoiMask = convhull(curSaSdT.mapIds(:, 1), curSaSdT.mapIds(:, 2));
    curRoiMask = curSaSdT.mapIds(curRoiMask, :);
    
    curRoiMask = poly2mask( ...
        curRoiMask(:, 1), curRoiMask(:, 2), ...
        curImSize(1), curImSize(2));
    
    % HACKHACKHACK(amotta): Generating the convex hull for a set of points
    % using `convhull` and then converting this hull to a mask using
    % `poly2mask` does not actually produce a mask that contains all of the
    % points from which it was derived. MATLAB is unbelievable...
    curRoiMask(curSaSdT.mapIdx) = true;
    
    curCtrlConfigs = struct('synIds', {}, 'title', {}, 'tag', {});
    curCtrlConfigs(1).synIds = curSaSdConfig.synIdPairs(:);
    curCtrlConfigs(1).title = 'SASD';
    
    %% Heat map
    for curCtrlConfig = curCtrlConfigs
        curKvPairs = { ...
            'xLim', curLimX, 'xAxis', curAxisX, ...
            'yLim', curLimY, 'yScale', curScaleY, ...
            'mapSize', curImSize, 'method', curMethod};
        
       [curSaSdMap, curBw] = ...
            connectEM.Consistency.densityMap( ...
                asiT.area(curSaSdConfig.synIdPairs), curKvPairs{:});
        
        %{
        % NOTE(amotta): Uncomment to simulate a null model in which random
        % synapse pairs are sampled from within the same connection type
        % (i.e., axon and target class pairs.
       [~, ~, curIds] = unique( ...
           curSaSdT(:, {'axonClass', 'targetClass'}), 'rows');
        curSaSdMap = ...
            connectEM.Consistency.nullDensityMaps( ...
                curSaSdT.areas, 'asiGroups', repmat(curIds, 1, 2), ...
                curKvPairs{:}, 'bandWidth', curBw, 'numMaps', 5000);
        curSaSdMap = mean(curSaSdMap, 3);
        %}
        
        curCtrlMaps = ...
            connectEM.Consistency.nullDensityMaps( ...
                asiT.area(curCtrlConfig.synIds), curKvPairs{:}, ...
                'bandWidth', curBw, 'numMaps', 5000);

        %% Prepare for figure
        curCtrlMap = mean(curCtrlMaps, 3);
        
        curMax = max(max(curSaSdMap(:)), max(curCtrlMap(:)));
        curDiffMap = curSaSdMap - curCtrlMap;
        curMaxDiff = max(abs(curDiffMap(:)));
        
        curPvalMap = 1 - mean(curCtrlMaps < curSaSdMap, 3);
        curPvalThreshs = sort(curConfig.twoDimPValThreshs, 'descend');
        
        curPvalImg = -log10(min( ...
            1 - mean(curCtrlMaps < curSaSdMap, 3), ...
            1 - mean(curCtrlMaps > curSaSdMap, 3)));
        
        % NOTE(amotta): Detect statistically significant regions. Drop tiny
        % regions (with less than 100 pixels or less than 1 % of SASD
        % connections), which are most likely caused by outliers.
        curRegionMask = curPvalMap < curPvalThreshs(1);
        curRegionMask = curRegionMask & curRoiMask;
        curRegionMask = bwlabel(curRegionMask);
        
        curRegionProps = regionprops( ...
            curRegionMask, {'Area', 'Centroid', 'BoundingBox'}); %#ok
        curKeepRegionIds = find([curRegionProps.Area] >= 100);
        
        curKeepRegionIds = curKeepRegionIds(arrayfun( ...
            @(id) sum(curSaSdMap(curRegionMask(:) == id)), ...
            curKeepRegionIds) > 0.01);
        
        curRegionProps = curRegionProps(curKeepRegionIds);
       [~, curRegionMask] = ismember(curRegionMask, curKeepRegionIds);
        curSaSdT.regionId = curRegionMask(curSaSdT.mapIdx);
        
        curConfigTitle = sprintf( ...
            '%s (n = %d pairs)', curConfig.title, ...
            size(curSaSdConfig.synIdPairs, 1));
        curCtrlTitle = sprintf( ...
            'vs. random pairs of %s (n = %d)', ...
            curCtrlConfig.title, floor(numel(curCtrlConfig.synIds) / 2));
        
        %% Figure
        curFig = figure();
        curFig.Color = 'white';
        curFig.Position(3:4) = [1060, 970];

        curAx = subplot(2, 2, 1);
        imagesc(curAx, curSaSdMap);
        caxis(curAx, [0, curMax]);
        colormap(curAx, jet(256));
        
        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
        
        curAx = subplot(2, 2, 2);
        imagesc(curAx, curCtrlMap);
        caxis(curAx, [0, curMax]);
        colormap(curAx, jet(256));

        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = {'0', sprintf('%.3g', curMax)};
        
        curAx = subplot(2, 2, 3);
        curPValAx = curAx;
        
        imagesc(curAx, curPvalImg);
        colormap(curAx, jet(256));
        
        curBar = colorbar('peer', curAx);
        curBar.Ticks = curBar.Limits;
        curBar.TickLabels = arrayfun( ...
            @(val) sprintf('%.3g', val), ...
            curBar.Limits, 'UniformOutput', false);
        curBar.Label.String = '-log10(p-value)';
        
        for curRegionId = 1:numel(curRegionProps)
            curPos = curRegionProps(curRegionId).Centroid;
            
            text(curAx, ...
                curPos(1), curPos(2), num2str(curRegionId), ...
                'Color', 'white', 'FontWeight', 'bold');
        end
        
        curAx = subplot(2, 2, 4);
        hold(curAx, 'on');
        
        imagesc(curAx, curDiffMap);
        caxis(curAx, [-1, +1] * curMaxDiff);
        colormap(curAx, jet(256));
        
        for curPvalThresh = curPvalThreshs
            contour(curAx, ...
                curRegionMask & ...
                curPvalMap < curPvalThresh, ...
                true, 'LineColor', 'black');
        end
        
        curBar = colorbar('peer', curAx);
        curBar.Ticks = [ ...
            curBar.Limits(1), ...
            mean(curBar.Limits), ...
            curBar.Limits(end)];
        curBar.TickLabels = { ...
            sprintf('%.3g', -curMaxDiff), '0', ...
            sprintf('%.3g', +curMaxDiff)};
        
        set( ...
            findobj(curFig.Children, 'Type', 'ColorBar'), ...
            'Location', 'EastOutside', ...
            'TickDirection', 'out', ...
            'Box', 'off');

        curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
        arrayfun(@(ax) hold(ax, 'on'), curAxes);
        
        for curMinPrctile = curMinPrctiles
            curAreaThresh = prctile( ...
                curSaSdT.areas(:), curMinPrctile);
            
            for curAx = curAxes
                contour(curAx, ...
                    curMinMap > curAreaThresh, true, ...
                    'LineStyle', '--', 'Color', 'white');
            end
        end

        curTickIdsX = 1 + floor((curImSize(2) - 1) * ...
            (curTicksX - curLimX(1)) / (curLimX(2) - curLimX(1)));
        curTickLabelsX = arrayfun( ...
            @num2str, curTicksX, 'UniformOutput', false);
        
        curTickIdsY = 1 + floor((curImSize(1) - 1) * ...
            (curTicksY - curLimY(1)) / (curLimY(2) - curLimY(1)));
        curTickLabelsY = arrayfun( ...
            @num2str, curTicksY, 'UniformOutput', false);

        set(curAxes, ...
            'Box', 'off', ...
            'TickDir', 'out', ...
            'YDir', 'normal', ...
            'YTick', curTickIdsY, ...
            'YTickLabels', curTickLabelsY, ...
            'YLim', [1, curImSize(2)], ...
            'XTick', curTickIdsX, ...
            'XTickLabels', curTickLabelsX, ...
            'XLim', [1, curImSize(1)], ...
            'PlotBoxAspectRatio', [1, 1, 1], ...
            'DataAspectRatioMode', 'auto');
        
        arrayfun(@(ax) xlabel(ax, curLabelX), curAxes);
        arrayfun(@(ax) ylabel(ax, curLabelY), curAxes);

        annotation( ...
            curFig, ...
            'textbox', [0, 0.9, 1, 0.1], ...
            'String', { ...
                info.filename; info.git_repos{1}.hash; ...
                curConfigTitle; curCtrlTitle}, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
        
        %% Evaluation
        curTitle = cell(numel(curRegionProps), 1);
        curFracsToPercsStr = @(f) strjoin(arrayfun( ...
            @(v) num2str(v, '%.1f%%'), 100 * f, ...
            'UniformOutput', false), ', ');
        
        for curRegionId = 1:numel(curRegionProps)
            curAreaRange = curRegionProps(curRegionId).BoundingBox;
            curAreaRange = curAreaRange(2) + [0, curAreaRange(4)];
            
            curAreaRange = ...
                curLimY(1) + diff(curLimY) ...
             .* curAreaRange / curImSize(1);
            curAreaRange = 10 .^ curAreaRange;
            
            curFracs = nan(numel(curPvalThreshs), 3);
            for curPvalIdx = 1:numel(curPvalThreshs)
                curPvalThresh = curPvalThreshs(curPvalIdx);
                
                curMask = curRegionMask == curRegionId;
                curMask = curMask & (curPvalMap < curPvalThresh);
                
                curFracs(curPvalIdx, 1) = sum(curSaSdMap(curMask));
                curFracs(curPvalIdx, 2) = sum(curCtrlMap(curMask));
                curFracs(curPvalIdx, 3) = mean(curMask(curSaSdT.mapIdx));
            end
            
            fprintf('Region %d\n', curRegionId);
            fprintf('* Average area: %.2f - %.2f µm²\n', curAreaRange);
            fprintf('* Upper bound (kde): %.1f - %.1f %%\n', 100 * curFracs(:, 1));
            fprintf('* Surplus (kde): %.1f - %.1f %%\n', 100 * diff(flip(curFracs(:, 1:2))));
            fprintf('* Upper bound (points): %.1f - %.1f %%\n', 100 * curFracs(:, 3));
            fprintf('\n');
            
            curTitle{curRegionId} = sprintf( ...
                'Region %d: %s', curRegionId, ...
                curFracsToPercsStr(curFracs(:, 3)));
        end
        
        curTitle = [{'Significance regions'}; curTitle]; %#ok
        title(curPValAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
        
        %% Evaluate pairs with small and large synapses
        %{
        curFig = figure();
       
        subplot(2, 3, 1);
        curAx = curFig.Children(1);
        
        curX = curMinMap(:);
       [curX, curIds] = sort(curX, 'ascend', 'MissingPlacement', 'last');
        
        hold(curAx, 'on');
        plot(gca, log10(curX), cumsum(curCtrlMap(curIds)));
        plot(gca, log10(curX), cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'log10(Max ASI area of small synapse [µm²])');
        xlim(curAx, curLimY);
        
        subplot(2, 3, 1 + 3);
        curAx = curFig.Children(1);
        
        % NOTE(amotta): Reuses `curX` from above
        curBinEdges = [-inf; sort(curSaSdT.areas(:)); +inf];
        curX = discretize(curX, curBinEdges);
        curX = curX / (numel(curBinEdges) - 1);
        
        hold(curAx, 'on');
        plot(gca, curX, cumsum(curCtrlMap(curIds)));
        plot(gca, curX, cumsum(curSaSdMap(curIds)));
        
        ylabel(curAx, 'Fraction of synapse pairs');
        xlabel(curAx, 'Max ASI area of small synapse [quantile]');
        xlim(curAx, [0, 1]);
        
        subplot(2, 3, 2);
        curAx = curFig.Children(1);
        
        curX = curMaxMap(:);
       [curX, curIds] = sort(curX, 'descend', 'MissingPlacement', 'last');
       
        hold(curAx, 'on');
        plot(curAx, log10(curX), cumsum(curCtrlMap(curIds)));
        plot(curAx, log10(curX), cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'log10(Min ASI area of large synapse [µm²])');
        xlim(curAx, curLimY);
        
        subplot(2, 3, 2 + 3);
        curAx = curFig.Children(1);
        
        % NOTE(amotta): Reuses `curX` from above
        curBinEdges = [-inf; sort(curSaSdT.areas(:)); +inf];
        curX = discretize(curX, curBinEdges);
        curX = curX / (numel(curBinEdges) - 1);
        
        hold(curAx, 'on');
        plot(gca, curX, cumsum(curCtrlMap(curIds)));
        plot(gca, curX, cumsum(curSaSdMap(curIds)));
        
        xlabel(curAx, 'Min ASI area of large synapse [quantile]');
        xlim(curAx, [0, 1]);
        
        subplot(2, 3, 3);
        curAx = curFig.Children(1);
        
        curX = linspace(curLimY(1), curLimY(2), curImSize(1));
       [curY, curX] = ndgrid(10 .^ curX, 10 .^ curX);
       
        curMap = arrayfun( ...
            @(x, y) sum(curSaSdMap( ...
                curMinMap(:) <= x ...
              & curMaxMap(:) >= y)), ...
            curX, curY);
        
        imagesc(curAx, curMap);
        caxis(curAx, [0, 1]);
        
        curPos = curAx.Position;
        curCbar = colorbar();
        curCbar.Label.String = 'Fraction of synapse pairs';
        curAx.Position = curPos;
        
        curTickIds = 1 + floor((curImSize(1) - 1) * ...
            (curTicksY - curLimY(1)) / (curLimY(2) - curLimY(1)));
        curTickLabels = arrayfun( ...
            @num2str, curTicksY, 'UniformOutput', false);
        
        xticks(curAx, curTickIds); xticklabels(curAx, curTickLabels);
        yticks(curAx, curTickIds); yticklabels(curAx, curTickLabels);
        
        xlabel(curAx, 'log10(Max ASI area of small synapse [µm²])');
        ylabel(curAx, 'log10(Min ASI area of large synapse [µm²])');
        
        curLines = findobj(curFig, 'type', 'line');
        set(curLines, 'LineWidth', 2);
        
        curAxes = flip(findobj(curFig, 'type', 'axes'));
        axis(curAxes, 'xy');
        axis(curAxes, 'square');
        
        curAx = curAxes(2);
        curPos = curAx.Position;
        curLeg = {'Random pairs', 'Same-axon same-dendrite pairs'};
        curLeg = legend(curAx, curLeg);
        curLeg.Location = 'SouthOutside';
        curAx.Position = curPos;
        
        annotation( ...
            curFig, 'textbox', [0, 0.9, 1, 0.1], ...
            'String', {curConfigTitle; curCtrlTitle});
        
        curFig.Position(3:4) = [1100, 770];
        connectEM.Figure.config(curFig, info);
        %}
        
        %% Evaluate sensitivity to p-value threshold
        %{
        curCtrlCounts = curCtrlMaps > curSaSdMap;
        curCtrlCounts = sum(curCtrlCounts, 3);
        
        curRegOff = 0;
        curRegMap = zeros(size(curCtrlCounts));
        curRegFracs = zeros(0, size(curCtrlMaps, 3));
        
        for curThresh = 1:size(curCtrlMaps, 3)
            curRegs = curCtrlCounts < curThresh;
            curRegs = regionprops(curRegs, 'PixelIdxList');
            curRegs = {curRegs.PixelIdxList};
            
            % Sort by size (for cosmetic reasons)
           [~, curIds] = sort(cellfun(@numel, curRegs));
            curRegs = curRegs(curIds);
            
            curRegIds = cellfun( ...
                @(ids) setdiff(curRegMap(ids), 0), ...
                curRegs, 'UniformOutput', false);
            
            curRegNewIds = find(cellfun(@isempty, curRegIds));
            curRegNewIds = reshape(curRegNewIds, 1, []);
            
            for curRegIdx = curRegNewIds
                curRegOff = curRegOff + 1;
                curRegIds{curRegIdx} = curRegOff;
                
                curRegMap(curRegs{curRegIdx}) = curRegOff;
                curRegFracs(curRegOff, :) = nan;
            end
            
            for curRegIdx = 1:numel(curRegs)
                curRegFrac = sum(curSaSdMap(curRegs{curRegIdx}));
                curRegFracs(curRegIds{curRegIdx}, curThresh) = curRegFrac;
            end
        end
        
        curFig = figure();
        curAx = axes(curFig); %#ok

        curX = linspace(0, 1, size(curCtrlMaps, 3));
        curY = transpose(flip(curRegFracs, 1));
        plot(curAx, curX, curY, 'LineWidth', 2);
        
        xlim(curAx, [0, 0.1]);
        xlabel(curAx, 'p-value threshold');
        ylabel(curAx, 'Fraction of SASD pairs in region');
        title(curAx, curConfigTitle);
        
        curFig.Position(3:4) = [330, 260];
        connectEM.Figure.config(curFig, info);
        %}
        
        %% Scatter plot
        curFig = figure();
        curAx = axes(curFig); %#ok
        hold(curAx, 'on');
        
        curData = curSaSdT;
        curData.pVal = curPvalMap(curData.mapIdx);
        curData = sortrows(curData, 'pVal', 'descend');
        
        if isempty(curRegionColors)
            curColors = get(groot, 'defaultAxesColorOrder');
            curColors = curColors(1, :);
        else
            curPValEdges = [0, sort(curPvalThreshs), 1];
            
            curData.colorIdx = discretize(curData.pVal, curPValEdges);
            curData.colorIdx = sub2ind( ...
               [numel(curPValEdges) - 1, 1 + numel(curRegionColors)], ...
                curData.colorIdx, 1 + curData.regionId);
            
            curColors = cellfun( ...
                @(colors) cat(1, colors, zeros(1, 3)), ...
                curRegionColors, 'UniformOutput', false);
            curColors = vertcat( ...
                zeros(1 + numel(curColors), 3), curColors{:});
            curColors = curColors(curData.colorIdx, :);
        end
        
        scatter(curAx, curData.x, curData.y, 1, curColors, '.');
        
        if ~isempty(curRegionContourProps)
            for curPvalThresh = curPvalThreshs
               [~, curCont] = contour(curAx, ...
                    curRegionMask & (curPvalMap < curPvalThresh), ...
                    true, curRegionContourProps{:});
                
                % NOTE(amotta): The X and Y coordinates of the contour are
                % in pixels relative to the p-value map. Let's convert this
                % to the `curLimX` and `curLimY` coordinate system.
                curCont.XData = ...
                    curLimX(1) + diff(curLimX) * ...
                    (curCont.XData - 1) / (curImSize(2) - 1);
                curCont.YData = ...
                    curLimY(1) + diff(curLimY) * ...
                    (curCont.YData - 1) / (curImSize(1) - 1);
            end
        end

        xlim(curAx, curLimX);
        xticks(curAx, curTicksX);
        xlabel(curAx, curLabelX);

        ylim(curAx, curLimY);
        yticks(curAx, curTicksY);
        ylabel(curAx, curLabelY);

        axis(curAx, 'square');
        connectEM.Figure.config(curFig, info);
        curFig.Position(3:4) = 160;
        
        %% Spatial distribution
        %{
        
        curVxSize = param.raw.voxelSize;
        
        curRegionNames = { ...
            'Non-significant pairs', ...
            'Small low-CV pairs', ...
            'Large low-CV pairs'};
        
        % TODO(amotta): Add legend
        error('TODO: Add legend');
        
        curColors = get(groot, 'defaultAxesColorOrder');
        curColors = cat(1, zeros(1, 3), curColors(1:2, :));
        
        curFig = figure();
        
        for curDimId = 1:3
            curDimName = char(int8('X') + (curDimId - 1));
            
            curLims = param.bbox .* curVxSize(:);
            curLims = curLims(curDimId, :) / 10E3;
            curLims = [floor(curLims(1)), ceil(curLims(2))] * 10;
            
            % Position
            curSynPos = asiT.pos(curSaSdT.synIds(:), curDimId);
            curSynPos = curSynPos .* curVxSize(curDimId) / 1E3;
            curSynPos = reshape(curSynPos, size(curSaSdT.synIds));
            
            % Distance
            curSynDist = diff(curSynPos, 1, 2);
            curSynDist = sqrt(curSynDist .* curSynDist);
            
            % Statistics
            for curRegionId = 1:2
                curCtrlMask = curSaSdT.regionId == 0;
                curThisMask = curSaSdT.regionId == curRegionId;
                
               [~, curPval] = kstest2( ...
                   curSynDist(curThisMask), ...
                   curSynDist(curCtrlMask));
            end
            
            % Plotting
            for curRegionId = 0:2
                curColor = curColors(curRegionId + 1, :);
                curSaSdMask = curSaSdT.regionId == curRegionId;
                
                % Position
                curAx = subplot(2, 3, curDimId);
                hold(curAx, 'on');
                
                histogram( ...
                    curAx, curSynPos(curSaSdMask, :), ...
                    'Normalization', 'probability', ...
                    'BinEdges', curLims(1):10:curLims(2), ...
                    'DisplayStyle', 'stairs', ...
                    'EdgeColor', curColor);
                
                xlabel(curAx, sprintf('%s position [µm]', curDimName));
                ylabel(curAx, 'Probability');
                
                % Distance
                curAx = subplot(2, 3, 3 + curDimId);
                hold(curAx, 'on');
                
                histogram( ...
                    curAx, curSynDist(curSaSdMask), ...
                    'Normalization', 'probability', ...
                    'BinEdges', 0:1:30, 'DisplayStyle', 'stairs', ...
                    'EdgeColor', curColor);
                
                xlabel(curAx, sprintf('%s distance [µm]', curDimName));
                ylabel(curAx, 'Probability');
            end
        end
        
        connectEM.Figure.config(curFig, info);
        %}
    end
end
