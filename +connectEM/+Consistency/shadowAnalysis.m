% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiRunId = '20190227T082543';

info = Util.runInfo();
Util.showRunInfo(info);

%% Load axon-spine interfaces
[curDir, curAsiFile] = fileparts(connFile);
curAsiFile = sprintf('%s__%s_asiT.mat', curAsiFile, asiRunId);
curAsiFile = fullfile(curDir, curAsiFile);

asiT = load(curAsiFile);
asiT = asiT.asiT;

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Prepare data
clear cur*;

plotConfigs = struct('synIds', {}, 'title', {}, 'tag', {});

axonClasses = { ...
    'Exc', {'Corticocortical', 'Thalamocortical'}; ...
    'CC',  {'Corticocortical'}; ...
    'TC',  {'Thalamocortical'}};
targetClasses = { ...
    'All', categories(asiT.targetClass); ...
    'PD',  'ProximalDendrite'; ...
    'AD',  'ApicalDendrite'; ...
    'OD',  'OtherDendrite'};

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
% NOTE(amotta): p-value threshold at which large and small low-CV regions
% are not yet merged. (I've only looked at range up to p = 10 %.)
curPvalThreshs = [0.5, 1, 2, 3, 4, 5] / 100;

% Excitatory axons
plotConfigs(1, 1).twoDimPValThreshs = 0.026;
plotConfigs(2, 1).twoDimPValThreshs = inf;
plotConfigs(3, 1).twoDimPValThreshs = 0.015;
plotConfigs(4, 1).twoDimPValThreshs = inf;

% Corticocortical axons
plotConfigs(1, 2).twoDimPValThreshs = 0.052;
plotConfigs(2, 2).twoDimPValThreshs = inf;
plotConfigs(3, 2).twoDimPValThreshs = inf;
plotConfigs(4, 2).twoDimPValThreshs = inf;

% Thalamocortical axons
plotConfigs(1, 3).twoDimPValThreshs = 0.087;
plotConfigs(2, 3).twoDimPValThreshs = inf;
plotConfigs(3, 3).twoDimPValThreshs = inf;
plotConfigs(4, 3).twoDimPValThreshs = inf;

for curIdx = 1:numel(plotConfigs)
    curThreshs = plotConfigs(curIdx).twoDimPValThreshs;
    
    % Only consider p-value thresholds before LTP and LTD regions merger
    curThreshs = curPvalThreshs(curPvalThreshs < curThreshs);
    
    % Only consider lowest and highest p-value threshold
    curThreshs = curThreshs(unique([1, numel(curThreshs)]));
    plotConfigs(curIdx).twoDimPValThreshs = curThreshs;
end

%%
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

    curSaSdT.mapIdx = [ ...
        discretize(curSaSdT.y, linspace( ...
            curLimY(1), curLimY(2), curImSize(1) + 1)), ...
        discretize(curSaSdT.x, linspace( ...
            curLimX(1), curLimX(2), curImSize(2) + 1))];
    curSaSdT.mapIdx = sub2ind( ...
        curImSize, curSaSdT.mapIdx(:, 1), curSaSdT.mapIdx(:, 2));
    
    curCtrlConfigs = struct('synIds', {}, 'title', {}, 'tag', {});
    curCtrlConfigs(1).synIds = curSaSdConfig.synIdPairs(:);
    curCtrlConfigs(1).title = 'SASD';
    
    %% Heat map        
    curKvPairs = { ...
        'xLim', curLimX, 'xAxis', curAxisX, ...
        'yLim', curLimY, 'yScale', curScaleY, ...
        'mapSize', curImSize, 'method', curMethod};
    
    curFig = figure();
    curFigs = curFig([]);
    delete(curFig);
    
    curCtrlConfigIdx = 0;
    while curCtrlConfigIdx < numel(curCtrlConfigs)
        curCtrlConfigIdx = curCtrlConfigIdx + 1;
        curCtrlConfig = curCtrlConfigs(curCtrlConfigIdx);
        
        curSaSdSynIdPairs = all(ismember( ...
            curSaSdConfig.synIdPairs, curCtrlConfig.synIds), 2);
        curSaSdSynIdPairs = curSaSdConfig.synIdPairs(curSaSdSynIdPairs, :);
        
       [curSaSdMap, curBw] = ...
            connectEM.Consistency.densityMap( ...
                asiT.area(curSaSdSynIdPairs), curKvPairs{:});
        
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
        curRegionMask = bwlabel(curRegionMask);
        
        curRegionProps = regionprops( ...
            curRegionMask, {'Area', 'Centroid'}); %#ok
        curKeepRegionIds = find([curRegionProps.Area] >= 100);
        
        curKeepRegionIds = curKeepRegionIds(arrayfun( ...
            @(id) sum(curSaSdMap(curRegionMask(:) == id)), ...
            curKeepRegionIds) > 0.01);
        
        curRegionProps = curRegionProps(curKeepRegionIds);
       [~, curRegionMask] = ismember(curRegionMask, curKeepRegionIds);
        curSaSdT.regionId = curRegionMask(curSaSdT.mapIdx);
        
        curConfigTitle = sprintf('%s (n = %d pairs)', ...
            curConfig.title, size(curSaSdSynIdPairs, 1));
        curCtrlTitle = sprintf( ...
            'vs. random pairs of %s (n = %d)', ...
            curCtrlConfig.title, floor(numel(curCtrlConfig.synIds) / 2));
        
        %% Figure
        curFig = figure();
        curFigs(end + 1) = curFig; %#ok
        
        curFig.Color = 'white';
        curFig.Position(3:4) = [1060, 970];

        curAx = subplot(2, 2, 1);
        imagesc(curAx, curSaSdMap);
        caxis(curAx, [0, curMax]);
        colormap(curAx, jet(256));
        colorbar('peer', curAx);
        
        curAx = subplot(2, 2, 2);
        imagesc(curAx, curCtrlMap);
        caxis(curAx, [0, curMax]);
        colormap(curAx, jet(256));
        colorbar('peer', curAx);
        
        curAx = subplot(2, 2, 3);
        curPValAx = curAx;
        
        imagesc(curAx, curPvalImg);
        colormap(curAx, jet(256));
        colorbar('peer', curAx);
        
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
        colorbar('peer', curAx);
        
        for curPvalThresh = curPvalThreshs
            contour(curAx, ...
                curRegionMask & ...
                curPvalMap < curPvalThresh, ...
                true, 'LineColor', 'black');
        end        
        
        set( ...
            findobj(curFig.Children, 'Type', 'ColorBar'), ...
            'Location', 'EastOutside', ...
            'TickDirection', 'out', ...
            'Box', 'off');

        curAxes = reshape(flip(findobj(curFig, 'Type', 'Axes')), 1, []);
        arrayfun(@(ax) hold(ax, 'on'), curAxes);

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
        for curRegionId = 1:numel(curRegionProps)
            curFracs = nan(numel(curPvalThreshs), 2);
            for curPvalIdx = 1:numel(curPvalThreshs)
                curPvalThresh = curPvalThreshs(curPvalIdx);
                
                curMask = curRegionMask == curRegionId;
                curMask = curMask & (curPvalMap < curPvalThresh);
                
                curFracs(curPvalIdx, 1) = sum(curSaSdMap(curMask));
                curFracs(curPvalIdx, 2) = sum(curCtrlMap(curMask));
            end
            
            curTitle{curRegionId} = sprintf( ...
                'Region %d: %s', curRegionId, strjoin(arrayfun( ...
                    @(v) num2str(v, '%.1f%%'), 100 * curFracs(:, 1), ...
                    'UniformOutput', false), ', '));
        end
        
        curTitle = [{'Significance regions'}; curTitle]; %#ok
        title(curPValAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
        
        %% Prepare new control conditions
        if ~isscalar(curCtrlConfigs); continue; end
        
        for curRegionId = 1:numel(curRegionProps)
            curSynIds = curSaSdT.regionId == curRegionId;
            curSynIds = curSaSdT.synIds(curSynIds, :);
            curSynIds = setdiff(curCtrlConfig(1).synIds, curSynIds);
            
            curTitle = sprintf( ...
                '%s without region %d', ...
                curCtrlConfig(1).title, curRegionId);
            
            curNewConfig = curCtrlConfigs([]);
            curNewConfig(1).synIds = curSynIds(:);
            curNewConfig(1).title = curTitle;
            curNewConfig(1).tag = [];
            
            curCtrlConfigs(end + 1) = curNewConfig; %#ok
        end
    end
    
    %% Fix color scales across figures
    curAxes = arrayfun( ...
        @(f) findobj(f, 'Type', 'Axes'), ...
        curFigs, 'UniformOutput', false);
    curAxes = cat(2, curAxes{:});
    
    curClims = get(curAxes, {'Clim'});
    curClims = reshape(curClims, size(curAxes));
    
    curClims = [ ...
        min(cellfun(@min, curClims), [], 2), ...
        max(cellfun(@max, curClims), [], 2)];
    
    curClims = num2cell(curClims, 2);
    curClims = repmat(curClims, 1, numel(curFigs));
    
    set(curAxes(:), {'CLim'}, curClims(:));
end
