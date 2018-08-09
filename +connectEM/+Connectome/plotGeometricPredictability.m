% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
availFile = '/tmpscratch/amotta/l4/2018-07-18-surface-availability-for-connectome-v7-partially-split/axon-availability_v2.mat';

minSynPre = 10;
maxRadius = 50;
maxAvail = 0.7;

% Valid prediction methods are
% * predictTargetClassAvailability: The predicted innervation of a target
%   class (in fractional synapses) is equal to the fraction of the
%   available postsynaptic membranes being of the given target class.
%
% * predictUsingLinearRegressionOnTargetClassAvailability: The predicted
%   innervation of a target class is the result of a linear regression
%   based on the fractional surface availability of the given target class.
%
% * predictUsingLinearRegressionOnAllTargetClassAvailabilities: The
%   predicted innervation of a target class is the result of a multivariate
%   linear regression based on the fractional surface availability of all
%   target classes.
% * predictUsingMultivariateMultinomialLogisticRegression: The predicted
%   innervation of a target class is the result of a multivariate
%   multinomial logistic regression based on the fractional surface
%   availability of all target classes.
predictionMethod = 'predictUsingMultivariateMultinomialLogisticRegression';

% Set to radius (in µm) to run forward model to generate fake connectome
% and calibrate the geometric predictability analysis.
fakeRadius = [];

targetClasses = { ...
    'Somata', 'SO';
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'Other'};

targetTags = reshape(targetClasses(:, 2), 1, []);
targetClasses = reshape(targetClasses(:, 1), 1, []);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[conn, axonClasses] = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, axonClasses, 'minSynPre', minSynPre);
axonClasses = axonClasses(1:2);
    
avail = load(availFile);

%% Prepare data
[classConn, classIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Fix order of target classes
[~, classIds] = ismember(targetClasses, classIds);
classConn = classConn(:, classIds);

% Determine relative availabilities of target classes
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

% Define target classes we want to predict
[axonClasses.predictClasses] = deal(targetClasses);
[axonClasses.correctForBinoVar] = deal(false);

% NOTE(amotta): As discussed with MH on 08.08.2018, let's not predict
% targets that are not innervated by excitatory axons.
axonClasses(1).title = sprintf( ...
    '%s (without SOM and AIS)', axonClasses(1).title);
axonClasses(1).predictClasses = setdiff( ...
    targetClasses, {'Somata', 'AxonInitialSegment'});

% Duplicate excitatory axons
axonClasses(end + 1) = axonClasses(1);
axonClasses(end).correctForBinoVar = true;
axonClasses(end).title = sprintf( ...
    '%s (uncorrected)', axonClasses(end).title);

% Inhibitory axons without correction
axonClasses(end + 1) = axonClasses(2);
axonClasses(end).correctForBinoVar = true;
axonClasses(end).title = sprintf( ...
    '%s (uncorrected)', axonClasses(end).title);

% Cosmetics
curColors = get(groot, 'defaultAxesColorOrder');
axonClasses(1).styles = {'Color', curColors(1, :), 'LineStyle', '-'};
axonClasses(2).styles = {'Color', curColors(2, :), 'LineStyle', '-'};
axonClasses(3).styles = {'Color', curColors(1, :), 'LineStyle', '--'};
axonClasses(4).styles = {'Color', curColors(2, :), 'LineStyle', '--'};
plotAxonClasses = 1:numel(axonClasses);

%% Build fake connectome for testing
if ~isempty(fakeRadius)
    fakeRadiusId = find(avail.dists == 1E3 * fakeRadius);
    fakeConn = availabilities(:, fakeRadiusId, :);
    fakeConn = transpose(squeeze(fakeConn));
    
    % Build fake connectome by multinomial sampling
    rng(0);
    fakeConn = cell2mat(cellfun( ...
        @(n, probs) mnrnd(n, probs), ...
        num2cell(sum(classConn, 2)), ...
        num2cell(fakeConn, 2), ...
        'UniformOutput', false));
    classConn = fakeConn;
    
    fakeAxonClassTitles = strcat( ...
        {'fake '}, {axonClasses.title}, ...
        sprintf(' (r_{fake} = %d µm)', fakeRadius));
   [axonClasses.title] = deal(fakeAxonClassTitles{:});
end

%% Plot mean availabilities / specificities
clear cur*;

curFig = figure();
curFig.Color = 'white';

curDists = avail.dists / 1E3;

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;

    curClassConn = classConn(curAxonIds, :);
    curSynCounts = sum(curClassConn, 2);

    for curTargetClassId = 1:numel(targetClasses)
        curAvails = availabilities(curTargetClassId, :, curAxonIds);
        curAvails = shiftdim(curAvails, 1);
        curMeanAvail = mean(curAvails, 2);
        
        curSpecs = curClassConn(:, curTargetClassId);
        curSpecs = curSpecs ./ sum(curClassConn, 2);
        curMeanSpec = mean(curSpecs);
        
        curAx = subplot( ...
            numel(axonClasses), numel(targetClasses), ...
            numel(targetClasses) * (curAxonClassId - 1) ...
          + curTargetClassId);
        hold(curAx, 'on');
        
        plot(curAx, curDists, curMeanAvail);
        plot(curAx, [0, maxRadius], repelem(curMeanSpec, 2));
        
        curAx.YLim(1) = 0;
        curAx.XLim = [0, maxRadius];
    end
end

curAx = flip(cat(1, curFig.Children));
set(curAx, 'TickDir', 'out');

for curIdx = 1:numel(targetClasses)
    title(curAx(curIdx), targetClasses{curIdx}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

for curIdx = 1:numel(axonClasses)
    ylabel( ...
        curAx(numel(targetClasses) * (curIdx - 1) + 1), ...
        {'Mean fraction of '; 'surface / synapses'});
end

xlabel(curAx(end - numel(targetClasses) + 1), 'Radius (µm)');
set(cat(1, curAx.Children), 'LineWidth', 2);

curAxPos = curAx(end).Position;
curLeg = legend(curAx(end), { ...
    'Availability', ...
    'Specificity'}, ...
    'Location', 'EastOutside');
curLeg.Box = 'off';
curAx(end).Position = curAxPos;

annotation(curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% Plot variances
clear cur*;

curFig = figure();
curFig.Color = 'white';

curDists = avail.dists / 1E3;

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;

    curClassConn = classConn(curAxonIds, :);
    curSynCounts = sum(curClassConn, 2);
    
    curTotalConnVar = 0;
    curTotalAvailBinoVar = zeros(size(curDists(:)));
    curTotalExplainedVar = zeros(size(curDists(:)));

    for curTargetClassId = 1:numel(targetClasses)
        curAvails = availabilities(curTargetClassId, :, curAxonIds);
        curAvails = shiftdim(curAvails, 1);

        curAvailBinoVar = curAvails .* (1 - curAvails);
        curAvailBinoVar = curAvailBinoVar ./ curSynCounts(:)';
        curAvailBinoVar = mean(curAvailBinoVar, 2);
        
        curSpecs = curClassConn(:, curTargetClassId);
        curSpecs = curSpecs ./ sum(curClassConn, 2);
        curConnVar = mean((curSpecs - mean(curSpecs)) .^ 2);
        
        curExplainedVar = nan(size(curDists));
        for curDistId = 1:numel(curDists)
            curPredIn = availabilities(:, curDistId, curAxonIds);
            curPredIn = transpose(squeeze(curPredIn));
            curPredIn(:, end + 1) = 1; %#ok
            
            curWarn = warning('off', 'MATLAB:rankDeficientMatrix');
            curPredParams = curPredIn \ curSpecs;
            warning(curWarn);
            
            curPreds = curPredIn * curPredParams;
            curPostPredVar = mean((curPreds - curSpecs) .^ 2);
            curDistExplainedVar = curConnVar - curPostPredVar;
            curExplainedVar(curDistId) = curDistExplainedVar;
        end
        
        curExplainableVar = curConnVar - curAvailBinoVar;
        
        % Accumulation
        curTotalConnVar = curTotalConnVar + curConnVar;
        curTotalAvailBinoVar = curTotalAvailBinoVar + curAvailBinoVar(:);
        curTotalExplainedVar = curTotalExplainedVar + curExplainedVar(:);
        
        curAx = subplot( ...
            numel(axonClasses), numel(targetClasses) + 1, ...
            (numel(targetClasses) + 1) * (curAxonClassId - 1) ...
          + curTargetClassId);
        hold(curAx, 'on');
        
        plot(curAx, curDists, curAvailBinoVar);
        plot(curAx, [0, maxRadius], repelem(curConnVar, 2));
        plot(curAx, curDists, curExplainableVar);
        plot(curAx, curDists, curExplainedVar);
    end
    
    curAxes = curFig.Children(1:numel(targetClasses));
    curMaxY = max(cat(2, curAxes.YLim));
    set(curAxes, 'YLim', [0, curMaxY]);
    
    % Plot total
    curAx = subplot( ...
        numel(axonClasses), numel(targetClasses) + 1, ...
        (numel(targetClasses) + 1) * (curAxonClassId - 1) ...
        + numel(targetClasses) + 1);
    hold(curAx, 'on');
    
    plot(curAx, curDists, curTotalAvailBinoVar);
    plot(curAx, [0, maxRadius], repelem(curTotalConnVar, 2));
    plot(curAx, curDists, curTotalConnVar(1) - curTotalAvailBinoVar);
    plot(curAx, curDists, curTotalExplainedVar);
end

curAx = flip(cat(1, curFig.Children));
set(curAx, 'XLim', [0, maxRadius], 'TickDir', 'out');

curTitles = cat(1, targetClasses(:), {'Total'});
for curIdx = 1:numel(curTitles)
    title(curAx(curIdx), curTitles{curIdx}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

for curIdx = 1:numel(axonClasses)
    ylabel( ...
        curAx((numel(targetClasses) + 1) * (curIdx - 1) + 1), ...
        'Variance');
end

xlabel(curAx(end - numel(targetClasses) + 1), 'Radius (µm)');
set(cat(1, curAx.Children), 'LineWidth', 2);

curAxPos = curAx(end).Position;
curLeg = legend(curAx(end), { ...
    'Binomial geometric variance', ...
    'Connectomic variance', ...
    'Explainable variance', ...
    'Explained variance'}, ...
    'Location', 'EastOutside');
curLeg.Box = 'off';
curAx(end).Position = curAxPos;

annotation(curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% Calculate variance introduced by multinomial model
clear cur*;

targetDistAxonMnVar = nan( ...
    numel(targetClasses), ...
    numel(avail.dists), ...
    numel(axonClasses));
targetAxonConnMnVar = nan( ...
    numel(targetClasses), ...
    numel(axonClasses));

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;

    curClassConn = classConn(curAxonIds, :);
    curSynCount = sum(curClassConn, 2);
    
    % Multinomial variability according to availabilities
    curAvails = availabilities(:, :, curAxonIds);
    curAvails = permute(curAvails, [3, 1, 2]);
    
    curAvailMnVar = curAvails .* (1 - curAvails) ./ curSynCount;
    curAvailMnVar = shiftdim(mean(curAvailMnVar, 1), 1);
    
    % Multinomial variability according to specificities
    curSpecs = curClassConn ./ curSynCount;
    curSpecMnVar = curSpecs .* (1 - curSpecs) ./ curSynCount;
    curSpecMnVar = shiftdim(mean(curSpecMnVar, 1), 1);
    
    targetDistAxonMnVar(:, :, curAxonClassId) = curAvailMnVar;
    targetAxonConnMnVar(:, curAxonClassId) = curSpecMnVar;
end

assert(~any(isnan(targetDistAxonMnVar(:))));
assert(~any(isnan(targetAxonConnMnVar(:))));

%% Calculate explainability
clear cur*;

curPredFunc = str2func(predictionMethod);
curWarnConf = warning('off', 'MATLAB:rankDeficientMatrix');

axonClassExplainability = nan( ...
    numel(avail.dists), numel(axonClasses));
axonTargetClassExplainability = nan( ...
    numel(targetClasses), numel(avail.dists), numel(axonClasses));

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curCorrect = curAxonClass.correctForBinoVar;
    curAxonIds = curAxonClass.axonIds;
    
    curPredictClasses = curAxonClass.predictClasses;
   [~, curPredictClassIds] = ismember(curPredictClasses, targetClasses);
    
    curConn = classConn(curAxonIds, :);
    curSynCounts = sum(curConn, 2);
    curConn = curConn ./ curSynCounts;
    curConn = curConn(:, curPredictClassIds);
    
    curVar = mean((curConn - mean(curConn, 1)) .^ 2, 1);
    
    tic();
    for curDistId = 1:numel(avail.dists)
        Util.progressBar(curDistId, numel(avail.dists));
        curAvails = availabilities(:, curDistId, curAxonIds);
        curAvails = transpose(squeeze(curAvails));
        curAvails(:, end + 1) = 1; %#ok
        
        curPred = curPredFunc(curConn, curAvails, curSynCounts);
        
        curVarUnexplainable = targetDistAxonMnVar( ...
            curPredictClassIds, curDistId, curAxonClassId);
        curVarUnexplainable = reshape(curVarUnexplainable, 1, []);
        curVarUnexplainable = curCorrect * curVarUnexplainable;
        
        % Per axon and target class
        curVarLeft = mean((curPred - curConn) .^ 2, 1);
        curVarLeft = max(curVarLeft - curVarUnexplainable, 0);
        
        curVarExplained = curVar - curVarLeft;
        curVarFracExplained = curVarExplained ./ curVar;
        
        axonTargetClassExplainability(curPredictClassIds, ...
            curDistId, curAxonClassId) = curVarFracExplained;
        
        % Per axon class
        curVarLeft = sum(curVarLeft);
        curVarExplained = sum(curVar) - curVarLeft;
        curVarFracExplained = curVarExplained / sum(curVar);
        
        axonClassExplainability( ...
            curDistId, curAxonClassId) = curVarFracExplained;
    end
end

warning(curWarnConf);
clear prevWarning;

%% Show availabilities for two axons
clear cur*;

% AD- and soma-specific axon, respectively
% TODO(amotta): Update to IDs in latest connectome!
curAxonIds = [];

if ~isempty(curAxonIds)
    curFig = figure();
    curFig.Color = 'white';

    curAx = axes(curFig);
    hold(curAx, 'on');
    axis(curAx, 'square');

    curAvail = availabilities(:, :, curAxonIds);
    curAvail = permute(curAvail, [2, 3, 1]);
    curAvail = reshape(curAvail, size(curAvail, 1), []);
    curAvail = transpose(curAvail);

    plot(curAx, avail.dists / 1E3, curAvail);

    curColors = curAx.ColorOrder(1:numel(targetClasses), :);
    curColors = num2cell(curColors, 2);

    curLines = flip(cat(1, curAx.Children));
    [curLines.LineWidth] = deal(2);
    [curLines(1:2:end).LineStyle] = deal('--');
    [curLines(1:2:end).Color] = deal(curColors{:});
    [curLines(2:2:end).Color] = deal(curColors{:});

    xlabel(curAx, 'r_{pred} (µm)');
    ylabel(curAx, 'Availability');

    curAx.XLim = [0, maxRadius];
    curAx.YLim = [0, maxAvail];
    curAx.TickDir = 'out';

    curLeg = legend( ...
        curLines(2:2:end), targetClasses, ...
        'Location', 'EastOutside');
    curLeg.Box = 'off';

    annotation(curFig, ...
        'textbox', [0, 0.9, 1, 0.1], ...
        'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
        'String', { ...
            info.filename; info.git_repos{1}.hash; ...
            'AD- (dashed) vs. soma-specific (solid) axon'});
end

%% Scatter plot and linear regression
clear cur*;

% Demonstrate how synapse fraction is predicted from availabilities
curSynCount = 10;
curRadius = 10;

curPad = 0.05;
curMinRange = 0.1;

curAxonIds = conn.axonMeta.synCount >= curSynCount;
curRadiusId = find(avail.dists == 1E3 * curRadius);

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [500, 310];

for curIdx = 1:numel(targetClasses)
    curSpecs = classConn(curAxonIds, :);
    curSpecs = curSpecs(:, curIdx) ./ sum(curSpecs, 2);
    curAvail = availabilities(curIdx, curRadiusId, curAxonIds);
    
    curRange = [min(curAvail), max(curAvail)];
    curRange(1) = floor(10 * curRange(1)) / 10;
    curRange(2) = ceil(10 * curRange(2)) / 10;
    
    curFit = fit(curAvail(:), curSpecs(:), 'poly1');
    curTitle = targetTags{curIdx};
    
    curAx = subplot(1, numel(targetClasses), curIdx);
    hold(curAx, 'on');
    
    scatter(curAx, curAvail(:), curSpecs(:), '.');
    plot(curAx, [0, 1], curFit([0, 1]), 'Color', 'black', 'LineWidth', 2);
    
    curAx.TickDir = 'out';
    curAx.XLim = curRange;
    curAx.XTick = curRange;
    
    curAx.Position(3) = 0.3 * diff(curRange);
    curAx.Position([2, 4]) = [0.2, 0.5];
    
    title(curAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
end

curAxes = flip(curFig.Children);
set(curAxes, 'YLim', [0, 1]);

curWidth = sum(arrayfun(@(a) a.Position(3), curAxes));
curWidth = curWidth + (numel(curAxes) - 1) * curPad;
curOffset = (1 - curWidth) / 2;

curAxes(1).Position(1) = curOffset;
for curIdx = 2:numel(curAxes)
    curAxes(curIdx).YAxis.Visible = 'off';
    curAxes(curIdx).Position(1) = ...
        sum(curAxes(curIdx - 1).Position([1, 3])) + curPad;
end

curAx = curAxes(1);
xlabel(curAx, 'Surface availability');
ylabel(curAx, 'Axonal synapse fraction');

curTitle = sprintf( ...
    'All axons with ≥ %d synapses. r_{pred} = %d µm', ...
    curSynCount, curRadius);

annotation(curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash; curTitle});

%% Plot R² per axon and dendrite class
clear cur*;

% Model: Linear combination of all availabilities
curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [1100, 550];

for curAxonClassIdx = 1:numel(plotAxonClasses)
    curAxonClassId = plotAxonClasses(curAxonClassIdx);
    curAxonClass = axonClasses(curAxonClassId);
    
    curPredictClasses = curAxonClass.predictClasses;
   [~, curPredictClassIds] = ismember(curPredictClasses, targetClasses);
    
    curAx = subplot(1, numel(plotAxonClasses), curAxonClassIdx);
    curAx.TickDir = 'out';
    
    axis(curAx, 'square');
    hold(curAx, 'on');
    
    for curTargetClassIdx = 1:numel(curPredictClassIds)
        curTargetClassId = curPredictClassIds(curTargetClassIdx);
        curColor = curAx.ColorOrder(curTargetClassId, :);
        
        curData = axonTargetClassExplainability( ...
            curTargetClassId, :, curAxonClassId);
        plot(curAx, ...
            avail.dists / 1E3, curData, ...
            'LineWidth', 2, 'Color', curColor);
    end
    
    title( ...
        curAx, curAxonClass.title, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

curAxes = cat(1, curFig.Children);
set(curAxes, ...
    'XLim', [0, maxRadius], ...
    'YLim', [0, 1]);
xlabel(curAxes(end), 'Radius (µm)');
ylabel(curAxes(end), 'R²');

legend(curFig.Children(1), ...
    arrayfun( ...
        @char, targetClasses, ...
        'UniformOutput', false), ...
    'Location', 'NorthEast', ...
    'Box', 'off');

annotation(curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', {predictionMethod; info.filename; info.git_repos{1}.hash});

%% Plot R² over classes
clear cur*;

% Plotting
curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [600, 590];

curAx = axes(curFig);
curAx.TickDir = 'out';
axis(curAx, 'square');
hold(curAx, 'on');

for curAxonClassId = plotAxonClasses
    curStyles = axonClasses(curAxonClassId).styles;
    curData = axonClassExplainability(:, curAxonClassId);
    plot(curAx, ...
        avail.dists / 1E3, curData, ...
        'LineWidth', 2, curStyles{:});
end

xlabel(curAx, 'Radius (µm)');
xlim(curAx, [0, maxRadius]);
ylabel(curAx, 'R²');
ylim(curAx, [0, 1]);

legend(curAx, ...
    {axonClasses(plotAxonClasses).title}, ...
    'Location', 'North', 'Box', 'off');

annotation(curFig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', {predictionMethod; info.filename; info.git_repos{1}.hash});
    
%% Geometric prediction models
function preds = predictTargetClassAvailability(~, avails, ~) %#ok
    preds = avails(:, 1:(end - 1));
end

function preds = predictUsingLinearRegressionOnTargetClassAvailability(conn, avails, ~) %#ok
    preds = nan(size(conn));
    for curIdx = 1:size(conn, 2)
        curAvails = avails(:, [curIdx, end]);
        curFit = curAvails \ conn(:, curIdx);
        preds(:, curIdx) = curAvails * curFit;
    end
end

function preds = predictUsingLinearRegressionOnAllTargetClassAvailabilities(conn, avails, ~) %#ok
    fit = avails \ conn;
    preds = avails * fit;
end

function preds = predictUsingMultivariateMultinomialLogisticRegression(conn, avails, synCounts) %#ok
    % HACK(amotta): Restore absolute synapse counts
    conn = round(conn .* synCounts);
    avails(:, end) = [];
    
    % HACK(amotta): Work around MATLAB throwing an error if multiple axons
    % have the same availabilities, but different synapses. This can be
    % solved by generating a single virtual axon with the pooled synapses.
   [X, ~, Y] = unique(avails, 'rows');
    Y = Y(:) + ((1:size(conn, 2)) - 1) * size(X, 1);
    Y = reshape(accumarray(Y(:), conn(:)), size(X, 1), []);
    
    warnConfig = warning('off');
    onCleanup(@() warning(warnConfig));
    
    fit = mnrfit(X, Y);
    preds = mnrval(fit, avails);
end
