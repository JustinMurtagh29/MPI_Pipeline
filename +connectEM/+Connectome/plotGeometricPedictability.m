% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear

%% Configuration
param = struct;
param.saveFolder = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/axon-avail-data.mat';
connName = 'connectome_axons_18_a_ax_spine_syn_clust';

minSynPre = 10;
info = Util.runInfo();

%% Loading data
conn = connectEM.Connectome.load(param, connName);
avail = load(availFile);

%% Prepare data
[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);

axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

synCounts = sum(classConn, 2);
axonsWithSynapses = find(synCounts);

axonClasses(1).axonIds = intersect( ...
    axonClasses(1).axonIds, axonsWithSynapses);
axonClasses(2).axonIds = intersect( ...
    axonClasses(2).axonIds, axonsWithSynapses);

% Determine relative availabilities of target classes
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

%% Calculate predictability
axonCount = size(classConn, 1);
distCount = numel(avail.dists);
classCount = numel(targetClasses);

bhattaDist = nan(axonCount, distCount);
distProbs = nan(axonCount, distCount);
distClassProbs = nan(axonCount, classCount, distCount);

for curDistIdx = 1:distCount
    curAvail = availabilities(:, curDistIdx, :);
    curAvail = transpose(squeeze(curAvail));
    
    % Calculate per-axon probabilities
    curProbs = mnpdf(classConn, curAvail);
    
    % Calculate per-synapse probabilities
    curProbs = curProbs .^ (1 ./ synCounts);
    distProbs(:, curDistIdx) = curProbs;
    
    % Calculate Bhattacharyya distance
    curBhattaDist = classConn ./ sum(classConn, 2);
    curBhattaDist = sqrt(curBhattaDist .* curAvail);
    curBhattaDist = -log(sum(curBhattaDist, 2));
    bhattaDist(:, curDistIdx) = curBhattaDist;
    
    % Calculate per-class probabilities
    for curClassIdx = 1:classCount
        curClassAvail = curAvail(:, curClassIdx);
        curClassAvail = [curClassAvail, (1 - curClassAvail)];
        
        curClassSyns = classConn(:, curClassIdx);
        curClassSyns = [curClassSyns, (synCounts - curClassSyns)];
        
        curClassProbs = mnpdf(curClassSyns, curClassAvail);
        curClassProbs = curClassProbs .^ (1 ./ synCounts);
        distClassProbs(:, curClassIdx, curDistIdx) = curClassProbs;
    end
end

%% Coefficient of determination
curAxonIds = axonClasses(2).axonIds;

curClassConn = classConn(curAxonIds, :);
curClassConn = curClassConn ./ sum(curClassConn, 2);

curSsTot = mean(curClassConn, 1);
curSsTot = sum((curClassConn - curSsTot) .^ 2, 1);

curSsRes = availabilities(:, :, curAxonIds);
curSsRes = permute(curSsRes, [3, 1, 2]);
curSsRes = sum((curClassConn - curSsRes) .^ 2, 1);
curSsRes = transpose(squeeze(curSsRes));

curRsq = 1 - (curSsRes ./ curSsTot);

% Plotting
fig = figure();
ax = axes(fig);
hold(ax, 'on');

for curIdx = 1:size(curRsq, 2)
    plot(ax, curRsq(:, curIdx), 'LineWidth', 2);
end

legend(arrayfun(@char, targetClasses, 'Uni', false));
% ylim([0, 0.5]);

%% Plotting
excProbs = distProbs(axonClasses(1).axonIds, :);
inhProbs = distProbs(axonClasses(2).axonIds, :);

fig = figure();
ax = axes(fig);
hold(ax, 'on');

plot(avail.dists, median(excProbs, 1), 'LineWidth', 2);
plot(avail.dists, median(inhProbs, 1), 'LineWidth', 2);

legend('Excitatory axons', 'Inhibitory axons');

ylim(ax, [0, 1]);

%% Plotting Bhattacharyya distance
excDists = bhattaDist(axonClasses(1).axonIds, :);
inhDists = bhattaDist(axonClasses(2).axonIds, :);

fig = figure();
ax = axes(fig);
hold(ax, 'on');

plot(avail.dists, median(excDists, 1), 'LineWidth', 2);
plot(avail.dists, median(inhDists, 1), 'LineWidth', 2);
ax.YLim(1) = 0;

%% Per-class plotting
fig = figure();

for curClassIdx = 1:classCount
    curAx = subplot(1, classCount, curClassIdx);
    hold(curAx, 'on');
    
    curExcProbs = axonClasses(1).axonIds;
    curExcProbs = distClassProbs(curExcProbs, curClassIdx, :);
    curExcProbs = median(squeeze(curExcProbs), 1);
    plot(curAx, avail.dists, curExcProbs, 'LineWidth', 2);
    
    curInhProbs = axonClasses(2).axonIds;
    curInhProbs = distClassProbs(curInhProbs, curClassIdx, :);
    curInhProbs = median(squeeze(curInhProbs), 1);
    plot(curAx, avail.dists, curInhProbs, 'LineWidth', 2);
    
    ylim(curAx, [0, 1]);
    
    title(curAx, ...
        char(targetClasses(curClassIdx)), ...
        'FontWeight', 'normal', 'FontSize', 10);
end