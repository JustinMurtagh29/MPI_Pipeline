% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
availFile = '/tmpscratch/amotta/l4/2018-04-27-surface-availability-connectome-v5/axon-availability_v2.mat';

minSynPre = 10;
maxRadius = 50;
info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
avail = load(availFile);

%% Prepare data
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

% Determine relative availabilities of target classes
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

% Count synapses per axon
% NOTE(amotta): Depending on the crazy stuff we're doing (e.g., excluding
% 'OtherDendrites') this value might be different from what's stored in
% conn.axonMeta.synCount.
synCounts = sum(classConn, 2);

%% Calculate 
% Approach: Let's use an axons observed specificities for a multinomial
% distribution and calculate its variance. Sum this up over all axons to
% get the expected sum of squares.
for curAxonClassId = 1:numel(axonClasses)
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    curAxonCount = numel(curAxonIds);
    
    curClassConn = classConn(curAxonIds, :);
    curSynCount = sum(curClassConn, 2);
    curClassConn = curClassConn ./ curSynCount;
    
    % The observed variance
    curVar = mean(curClassConn, 1);
    curVar = (curClassConn - curVar) .^ 2;
    curVar = sum(curVar(:)) / curAxonCount;
    
    % The variance of output variable i in a multinomail distribution is
    % n * p * (1 - p). The variance of the specificity (i.e., after
    % normalization) is thus p * (1 - p) / n.
    curMnVar = curClassConn .* (1 - curClassConn) ./ curSynCount;
    curMnVar = sum(curMnVar(:)) / curAxonCount;
    
    % Show result
    fprintf('%s\n', curAxonClass.title);
    fprintf('* Variance: %g\n', curVar);
    fprintf('* Multinomial variance: %g\n', curMnVar);
    fprintf('* Maximum explainable: %.2f %%\n', ...
        100 * (1 - curMnVar / curVar));
    fprintf('\n');
end

%% Plot R² per axon and dendrite class
% Model: Availability is synapse fraction
plotAxonClasses = 1:2;
plotTargetClasses = targetClasses;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1100, 550];

for curAxonClassIdx = 1:numel(plotAxonClasses)
    curAxonClassId = plotAxonClasses(curAxonClassIdx);
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    
    % Fractional class connectome
    curClassConn = classConn(curAxonIds, :);
    curClassConn = curClassConn ./ sum(curClassConn, 2);
    
    % Variance
    curSsTot = mean(curClassConn, 1);
    curSsTot = sum((curClassConn - curSsTot) .^ 2, 1);
    
    % Residuals
    curSsRes = availabilities(:, :, curAxonIds);
    curSsRes = permute(curSsRes, [3, 1, 2]);
    curSsRes = sum((curClassConn - curSsRes) .^ 2, 1);
    curSsRes = transpose(squeeze(curSsRes));
    
    % R²
    curRsq = 1 - (curSsRes ./ curSsTot);
    
    % Plotting
    curAx = subplot(1, numel(plotAxonClasses), curAxonClassIdx);
    curAx.TickDir = 'out';
    axis(curAx, 'square');
    hold(curAx, 'on');
    
    for curTargetClassIdx = 1:numel(plotTargetClasses)
        plot( ...
            curAx, avail.dists / 1E3, ...
            curRsq(:, curTargetClassIdx), ...
            'LineWidth', 2);
    end
    
    % Hint zero line
    plot(curAx, avail.dists([1, end]) / 1E3, [0, 0], 'k--');
    
    title( ...
        curAx, curAxonClass.title, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

[fig.Children.XLim] = deal([0, maxRadius]);
[fig.Children.YLim] = deal([-0.3, 0.5]);
xlabel(fig.Children(end), 'Radius (µm)');
ylabel(fig.Children(end), 'R²');

legend(fig.Children(1), ...
    arrayfun( ...
        @char, plotTargetClasses, ...
        'UniformOutput', false), ...
	'Location', 'NorthEast', ...
    'Box', 'off');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        'Prediction is availability'; ...
        info.filename; info.git_repos{1}.hash});

%% Plot R² per axon class
% Model: Availability is synapse fraction
plotAxonClasses = 1:2;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [600, 590];

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

for curAxonClassIdx = 1:numel(plotAxonClasses)
    curAxonClassId = plotAxonClasses(curAxonClassIdx);
    curAxonClass = axonClasses(curAxonClassId);
    curAxonIds = curAxonClass.axonIds;
    
    % Fractional class connectome
    curClassConn = classConn(curAxonIds, :);
    curClassConn = curClassConn ./ sum(curClassConn, 2);
    
    % Variance
    curSsTot = mean(curClassConn, 1);
    curSsTot = sum(sum((curClassConn - curSsTot) .^ 2));
    
    % Residuals
    curSsRes = availabilities(:, :, curAxonIds);
    curSsRes = permute(curSsRes, [3, 1, 2]);
    curSsRes = (curClassConn - curSsRes) .^ 2;
    curSsRes = reshape(sum(sum(curSsRes, 1), 2), [], 1);
    
    % R²
    curRsq = 1 - (curSsRes ./ curSsTot);
    
    % Plotting
    plot( ...
        ax, avail.dists / 1E3, ...
        curRsq, 'LineWidth', 2);
end

plot(ax, ...
    avail.dists([1, end]) / 1E3, [0, 0], ...
    'Color', 'black', 'LineStyle', '--');
legend(ax, ...
    {axonClasses(plotAxonClasses).title}, ...
    'Location', 'North', 'Box', 'off');

xlabel(ax, 'Radius (µm)');
xlim(ax, [0, maxRadius]);
ylabel(ax, 'R²');
ylim(ax, [-0.3, 0.5]);

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        'Prediction is availability'; ...
        info.filename; info.git_repos{1}.hash});

%% Plot R² per axon and dendrite class
% Model: Linear combination of all availabilities
plotAxonClasses = 1:2;

% Prepare output
rSq = nan( ...
    numel(targetClasses), ...
    numel(avail.dists), ...
    numel(plotAxonClasses));

% Calculate all R² values
for curAxonClassIdx = 1:numel(plotAxonClasses)
    curAxonClassId = plotAxonClasses(curAxonClassIdx);
    curAxonIds = axonClasses(curAxonClassId).axonIds;

    % Fractional connectome
    curClassConn = classConn(curAxonIds, :);
    curClassConn = curClassConn ./ sum(curClassConn, 2);

    curSsTot = mean(curClassConn, 1);
    curSsTot = sum((curClassConn - curSsTot) .^ 2, 1);

    for curDistIdx = 1:numel(avail.dists)
        curAvail = availabilities(:, curDistIdx, curAxonIds);
        curAvail = transpose(squeeze(curAvail));
        curAvail(:, end + 1) = 1; %#ok
        
        curCoefs = curAvail \ curClassConn;
        curPreds = curAvail * curCoefs;

        % Calculate sum of squares of residuals
        curSsRes = sum((curPreds - curClassConn) .^ 2, 1);
        curRSq = 1 - (curSsRes ./ curSsTot);
        
        rSq(:, curDistIdx, curAxonClassIdx) = curRSq;
    end
end

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
    
    for curTargetClassIdx = 1:numel(plotTargetClasses)
        plot( ...
            curAx, avail.dists / 1E3, ...
            rSq(curTargetClassIdx, :, curAxonClassIdx), ...
            'LineWidth', 2);
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
        @char, plotTargetClasses, ...
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

% Prepare output
rSq = nan( ...
    numel(avail.dists), ...
    numel(plotAxonClasses));

% Calculate all R² values
for curAxonClassIdx = 1:numel(plotAxonClasses)
    curAxonClassId = plotAxonClasses(curAxonClassIdx);
    curAxonIds = axonClasses(curAxonClassId).axonIds;

    % Fractional connectome
    curClassConn = classConn(curAxonIds, :);
    curClassConn = curClassConn ./ sum(curClassConn, 2);

    curSsTot = mean(curClassConn, 1);
    curSsTot = sum(sum((curClassConn - curSsTot) .^ 2));

    for curDistIdx = 1:numel(avail.dists)
        curAvail = availabilities(:, curDistIdx, curAxonIds);
        curAvail = transpose(squeeze(curAvail));
        curAvail(:, end + 1) = 1; %#ok
        
        curCoefs = curAvail \ curClassConn;
        curPreds = curAvail * curCoefs;

        % Calculate sum of squares of residuals
        curSsRes = sum((curPreds(:) - curClassConn(:)) .^ 2);
        curRSq = 1 - (curSsRes ./ curSsTot);
        
        rSq(curDistIdx, curAxonClassIdx) = curRSq;
    end
end

% Plotting
fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [600, 590];

ax = axes(fig);
ax.TickDir = 'out';
axis(ax, 'square');
hold(ax, 'on');

plot(ax, avail.dists / 1E3, rSq, 'LineWidth', 2);

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
