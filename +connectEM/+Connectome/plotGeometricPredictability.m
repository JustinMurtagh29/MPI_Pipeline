% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
availFile = '/tmpscratch/amotta/l4/2018-04-27-surface-availability-connectome-v5/axon-availability_v2.mat';

maxRadius = 50;
maxAvail = 0.7;
plotAxonClasses = 1:2;

% Set to radius (in µm) to run forward model to generate fake connectome
% and calibrate the geometric predictability analysis.
fakeRadius = [];

targetClasses = { ...
    'Somata', ...
    'ProximalDendrite', ...
    'SmoothDendrite', ...
    'ApicalDendrite', ...
    'AxonInitialSegment', ...
    'OtherDendrite'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = connectEM.Connectome.load(param, connFile);
avail = load(availFile);

%% Prepare data
% Rename whole cells to proximal dendrites
% TODO(amotta): Use `connectEM.Connectome.prepareForSpecificityAnalysis`
% once availability calculation was updated accordingly!
wcMask = conn.denMeta.targetClass == 'WholeCell';
conn.denMeta.targetClass(wcMask) = 'ProximalDendrite';

wcMask = avail.targetClasses == 'WholeCell';
avail.targetClasses(wcMask) = 'ProximalDendrite';

% Add axon class tags
axonClasses(1).tag = 'Exc';
axonClasses(2).tag = 'Inh';
axonClasses(3).tag = 'TC';
axonClasses(4).tag = 'CC';

[classConn, classIds] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Fix order of target classes
[~, classIds] = ismember(targetClasses, classIds);
classConn = classConn(:, classIds);

% Determine relative availabilities of target classes
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

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

%% Calculate amount of variance introduced by multinomial model
targetDistAxonMnVar = nan( ...
    numel(targetClasses), ...
    numel(avail.dists), ...
    numel(axonClasses));

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;

    curSynCount = classConn(curAxonIds, :);
    curSynCount = sum(curSynCount, 2);
    
    curAvails = availabilities(:, :, curAxonIds);
    curAvails = permute(curAvails, [3, 1, 2]);
    
    curMnVar = curAvails .* (1 - curAvails) ./ curSynCount;
    curMnVar = shiftdim(mean(curMnVar, 1), 1);
    
    targetDistAxonMnVar(:, :, curAxonClassId) = curMnVar;
end

assert(~any(isnan(targetDistAxonMnVar(:))));

%% Testing
%{
fakeRadius = 20;

fakeRadiusId = find(avail.dists == 1E3 * fakeRadius);
fakeConn = availabilities(:, fakeRadiusId, :);
fakeConn = permute(fakeConn, [3, 1, 2]);

% Build fake connectome by multinomial sampling
rng(0);
fakeConn = cell2mat(cellfun( ...
    @(n, p) mnrnd(n, p), ...
    num2cell(sum(classConn, 2)), ...
    num2cell(fakeConn, 2), ...
    'UniformOutput', false));
fakeConn = fakeConn ./ sum(fakeConn, 2);

%% Single-target test case
curAxonClassId = 1;
curAxonClass = axonClasses(curAxonClassId);
curAxonIds = curAxonClass.axonIds;

curConn = fakeConn(curAxonIds, :);
curVar = (curConn - mean(curConn, 2)) .^ 2;
curVar = mean(curVar(:));

curVarFracsExplained = nan(size(avail.dists));
for curDistId = 1:numel(avail.dists)
    curAvails = availabilities(:, curDistId, curAxonIds);
    curAvails = transpose(squeeze(curAvails));
    curAvails(:, end + 1) = 1; %#ok
    
    curFit = curAvails \ curConn;
    curPred = curAvails * curFit;
    
    curVarLeft = mean(sum((curPred - curConn) .^ 2, 2), 1);
    curVarExplained = curVar - curVarLeft;
    
    curVarUnexplainable = sum( ...
        targetDistAxonMnVar(:, curDistId, curAxonClassId));
    curVarExplainable = curVar - curVarUnexplainable;
    
    curVarFracExplained = curVarExplained / curVarExplainable;
    curVarFracsExplained(curDistId) = curVarFracExplained;
end

fig = figure;
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');

plot(avail.dists / 1E3, curVarFracsExplained);

legend('Actual variance', 'Binomial variance');
ax.YLim(1) = 0;
%}

%% Show availabilities for two axons
% AD- and soma-specific axon, respectively
curAxonIds = [23314, 628];

fig = figure();
fig.Color = 'white';

ax = axes(fig);
hold(ax, 'on');
axis(ax, 'square');

curAvail = availabilities(:, :, curAxonIds);
curAvail = permute(curAvail, [2, 3, 1]);
curAvail = reshape(curAvail, size(curAvail, 1), []);
curAvail = transpose(curAvail);

plot(ax, avail.dists / 1E3, curAvail);

colors = ax.ColorOrder(1:numel(targetClasses), :);
colors = num2cell(colors, 2);

lines = flip([ax.Children]);
[lines.LineWidth] = deal(2);
[lines(1:2:end).LineStyle] = deal('--');
[lines(1:2:end).Color] = deal(colors{:});
[lines(2:2:end).Color] = deal(colors{:});

xlabel(ax, 'r_{pred} (µm)');
ylabel(ax, 'Availability');

ax.XLim = [0, maxRadius];
ax.YLim = [0, maxAvail];
ax.TickDir = 'out';

leg = legend( ...
    lines(2:2:end), targetClasses, ...
    'Location', 'EastOutside');
leg.Box = 'off';

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', { ...
        info.filename; info.git_repos{1}.hash; ...
        'AD- (dashed) vs. soma-specific (solid) axon'});

%% Scatter plot and linear regression
% Demonstrate how synapse fraction is predicted from availabilities
curSynCount = 10;
curRadius = 10;

curAxonIds = conn.axonMeta.synCount >= curSynCount;
curRadiusId = find(avail.dists == 1E3 * curRadius);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1150, 350];

for curIdx = 1:numel(targetClasses)
    curSpecs = classConn(curAxonIds, :);
    curSpecs = curSpecs(:, curIdx) ./ sum(curSpecs, 2);
    curAvail = availabilities(curIdx, curRadiusId, curAxonIds);
    
   [curFit, curGof] = fit(curAvail(:), curSpecs(:), 'poly1');
   
    curTitle = { ...
        char(targetClasses(curIdx)); ...
        sprintf('r² = %g', curGof.rsquare)};
    
    curAx = subplot(1, numel(targetClasses), curIdx);
    axis(curAx, 'square');
    hold(curAx, 'on');
    
    scatter(curAx, curAvail(:), curSpecs(:), '.');
    plot(curAx, [0, 1], curFit([0, 1]), 'Color', 'black');
    
    curAx.XLim = [0, 1];
    curAx.YLim = [0, 1];
    curAx.TickDir = 'out';
    title(curAx, curTitle, 'FontWeight', 'normal', 'FontSize', 10);
end

curAx = fig.Children(end);
xlabel(curAx, 'Availability');
ylabel(curAx, 'Synapse fraction');

curTitle = sprintf( ...
    'All axons with ≥ %d synapses. r_{pred} = %d µm', ...
    curSynCount, curRadius);

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash; curTitle});

%% Plot R² per axon and dendrite class
% Model: Linear combination of all availabilities
axonClassExplainability = nan( ...
    numel(avail.dists), numel(axonClasses));
axonTargetClassExplainability = nan( ...
    numel(targetClasses), numel(avail.dists), numel(axonClasses));

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    
    curConn = classConn(curAxonIds, :);
    curConn = curConn ./ sum(curConn, 2);
    
    curVar = mean((curConn - mean(curConn, 1)) .^ 2, 1);
    
    for curDistId = 1:numel(avail.dists)
        curAvails = availabilities(:, curDistId, curAxonIds);
        curAvails = transpose(squeeze(curAvails));
        curAvails(:, end + 1) = 1; %#ok
        
        curFit = curAvails \ curConn;
        curPred = curAvails * curFit;
        
        
        % Per axon and target class
        curVarLeft = mean((curPred - curConn) .^ 2, 1);
        curVarExplained = curVar - curVarLeft;
        
        curVarUnexplainable = ...
            targetDistAxonMnVar(:, curDistId, curAxonClassId);
        curVarUnexplainable = reshape(curVarUnexplainable, 1, []);
        
        curVarExplainable = curVar - curVarUnexplainable;
        curVarFracExplained = curVarExplained ./ curVarExplainable;
        
        axonTargetClassExplainability(:, ...
            curDistId, curAxonClassId) = curVarFracExplained;
        
        
        % Per axon class
        curVarLeft = sum(curVarLeft);
        curVarExplained = sum(curVar) - curVarLeft;
        
        curVarUnexplainable = sum(curVarUnexplainable);
        curVarExplainable = sum(curVar) - curVarUnexplainable;
        curVarFracExplained = curVarExplained / curVarExplainable;
        
        axonClassExplainability( ...
            curDistId, curAxonClassId) = curVarFracExplained;
    end
end

%%
% Plotting
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1100, 550];

for curAxonClassIdx = 1:numel(plotAxonClasses)
    curAxonClassId = plotAxonClasses(curAxonClassIdx);
    curAxonClass = axonClasses(curAxonClassId);
    
    curAx = subplot(1, numel(plotAxonClasses), curAxonClassIdx);
    curAx.TickDir = 'out';
    
    axis(curAx, 'square');
    hold(curAx, 'on');
    
    for curTargetClassIdx = 1:numel(targetClasses)
        curData = axonTargetClassExplainability( ...
            curTargetClassIdx, :, curAxonClassIdx);
        plot(curAx, avail.dists / 1E3, curData, 'LineWidth', 2);
    end
    
    title( ...
        curAx, curAxonClass.title, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

[fig.Children.XLim] = deal([0, maxRadius]);
[fig.Children.YLim] = deal([0, 1]);
xlabel(fig.Children(end), 'Radius (µm)');
ylabel(fig.Children(end), 'R²');

legend(fig.Children(1), ...
    arrayfun( ...
        @char, targetClasses, ...
        'UniformOutput', false), ...
    'Location', 'NorthEast', ...
    'Box', 'off');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        'Prediction by linear combination of availabilities'; ...
        info.filename; info.git_repos{1}.hash});

%% Plot R² over all classes
plotAxonClasses = 1:2;

% Plotting
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [600, 590];

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

for curAxonClassId = plotAxonClasses
    curData = axonClassExplainability(:, curAxonClassId);
    plot(ax, avail.dists / 1E3, curData, 'LineWidth', 2);
end

xlabel(ax, 'Radius (µm)');
xlim(ax, [0, maxRadius]);
ylabel(ax, 'R²');
ylim(ax, [0, 1]);

legend(ax, ...
    {axonClasses(plotAxonClasses).title}, ...
    'Location', 'North', 'Box', 'off');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        'Prediction by linear combination of availabilities'; ...
        info.filename; info.git_repos{1}.hash});

%% Plot synapse data
axonTargetClassMean = nan(numel(axonClasses), numel(targetClasses)); 
axonTargetClassVar = nan(numel(axonClasses), numel(targetClasses));

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    
    curSynCount = classConn(curAxonIds, :);
    curSynCount = curSynCount ./ sum(curSynCount, 2);
    
    curMean = mean(curSynCount, 1);
    curVar = (curSynCount - curMean) .^ 2;
    curVar = sum(curVar, 1) ./ numel(curAxonIds);
    
    axonTargetClassMean(curAxonClassId, :) = curMean;
    axonTargetClassVar(curAxonClassId, :) = curVar;
end

legends = arrayfun(@char, targetClasses, 'UniformOutput', false);

fig = figure();
fig.Color = 'white';

ax = subplot(3, 1, 1);
bar(ax, axonTargetClassMean, 'stacked');

ax.TickDir = 'out';
legend(ax, legends, 'Location', 'EastOutside');
xticklabels(ax, {axonClasses.tag});
ylabel(ax, 'Synapse fraction');

ax = subplot(3, 1, 2);
bar(ax, axonTargetClassVar, 'stacked');

ax.TickDir = 'out';
legend(ax, legends, 'Location', 'EastOutside');
xticklabels(ax, {axonClasses.tag});
ylabel(ax, 'Variance');

ax = subplot(3, 1, 3);
bar(ax, axonTargetClassBinoVar, 'stacked');

ax.TickDir = 'out';
legend(ax, legends, 'Location', 'EastOutside');
xticklabels(ax, {axonClasses.tag});
ylabel(ax, 'Binomial variance');

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Plot availability data
curLimY = 0.25;
curVarGain = 1;
curDists = avail.dists(:)' / 1E3;
curLegends = arrayfun( ...
    @char, targetClasses, ...
    'UniformOutput', false);

fig = figure();
fig.Color = 'white';

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    
    curAvails = availabilities(:, :, curAxonIds);
    curMean = mean(curAvails, 3);
    curVar = (curAvails - curMean) .^ 2;
    curVar = sum(curVar, 3) / numel(curAxonIds);
    
    curAx = subplot(numel(axonClasses), 1, curAxonClassId);
    hold(curAx, 'on');
    
    curMasterPlot = plot( ...
        curAx, curDists, curMean, 'LineWidth', 2);
    curPlot = plot(curAx, curDists, ...
        curMean + curVarGain * curVar, 'LineWidth', 0.5);
   [curPlot.Color] = deal(curMasterPlot.Color);
    curPlot = plot(curAx, curDists, ...
        curMean - curVarGain * curVar, 'LineWidth', 0.5);
   [curPlot.Color] = deal(curMasterPlot.Color);
   
   
    legend(curMasterPlot, curLegends, 'Location', 'EastOutside');
    title(curAx, curAxonClass.title, 'FontWeight', 'normal', 'FontSize', 10);
    
    curAx.TickDir = 'out';
    xlim(curAx, [0, maxRadius]);
    ylim(curAx, [0, curLimY]);
end

xlabel(curAx, 'Radius (µm)');
ylabel(curAx, 'Availabilities (mean ± std)');

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Synapse fraction vs. binomial variance scatter plot
curFigSize = [numel(axonClasses), numel(targetClasses)];

fig = figure();
fig.Color = 'white';

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    
    curSynCount = classConn(curAxonIds, :);
    curSynCount = sum(curSynCount, 2);
    curSynCount = curSynCount ./ curSynCount;
    
    for curTargetClassId = 1:numel(targetClasses)
        curProb = curSynCount(:, curTargetClassId);
        curVar = curProb .* (1 - curProb) ./ curSynCount;
        
        curAx = subplot( ...
            curFigSize(1), curFigSize(2), ...
            (curAxonClassId - 1) * curFigSize(2) + curTargetClassId);
        
        hold(curAx, 'on');
        scatter(curAx, curProb, curVar, '.');
        plot(curAx, [0, 1], repelem(mean(curVar), 1, 2));
        
        if curAxonClassId == 1
            title( ...
                curAx, char(targetClasses(curTargetClassId)), ...
                'FontWeight', 'normal', 'FontSize', 10);
        end
        
        axis(curAx, 'square');
    end
end

[fig.Children.YLim] = deal([0, 0.03]);
[fig.Children.XLim] = deal([0, 1]);

curAx = subplot( ...
    curFigSize(1), curFigSize(2), ...
    prod(curFigSize) - (curFigSize(2) - 1));
xlabel(curAx, 'Synapse probability');
ylabel(curAx, 'Binomial variance');

annotation( ...
    fig, 'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
