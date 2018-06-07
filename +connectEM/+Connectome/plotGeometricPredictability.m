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

%% Poking at geometric predictability model
fig = figure();
fig.Color = 'white';

curDists = avail.dists / 1E3;

for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;

    curClassConn = classConn(curAxonIds, :);
    curSynCounts = sum(curClassConn, 2);

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
        
        curAx = subplot( ...
            numel(axonClasses), numel(targetClasses), ...
            numel(targetClasses) * (curAxonClassId - 1) ...
          + curTargetClassId);
        hold(curAx, 'on');
        
        plot(curAx, curDists, curAvailBinoVar);
        plot(curAx, [0, maxRadius], repelem(curConnVar, 2));
        plot(curAx, curDists, curExplainableVar);
        plot(curAx, curDists, curExplainedVar);
        
        curAx.YLim(1) = 0;
        curAx.XLim = [0, maxRadius];
    end
end

ax = flip(cat(1, fig.Children));
set(ax, 'TickDir', 'out');

for curIdx = 1:numel(targetClasses)
    title(ax(curIdx), targetClasses{curIdx}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

plots = cat(1, ax.Children);
set(plots, 'LineWidth', 2);

axPos = ax(end).Position;
leg = legend(ax(end), { ...
    'Binomial geometric variance', ...
    'Connectomic variance', ...
    'Explainable variance', ...
    'Explained variance'}, ...
    'Location', 'EastOutside');
leg.Box = 'off';
ax(end).Position = axPos;

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% Calculate variance introduced by multinomial model
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
prevWarning = warning('off', 'MATLAB:rankDeficientMatrix');

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

warning(prevWarning);
clear prevWarning;

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
[fig.Children.YLim] = deal([-0.2, 1.2]);
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

%% Plot R² over classes
% Model: Linear combination of all availabilities
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
ylim(ax, [-0.2, 1.2]);

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
