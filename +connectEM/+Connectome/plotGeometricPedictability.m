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
synCounts = conn.axonMeta.synCount;

[classConn, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);
axonClasses = ...
    connectEM.Connectome.buildAxonClasses(conn, 'minSynPre', minSynPre);

% Determine relative availabilities of target classes
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

%% Calculate predictability
distCount = numel(avail.dists);
distProbs = nan(size(classConn, 1), distCount);

for curDistIdx = 1:distCount
    curAvail = availabilities(:, curDistIdx, :);
    curAvail = transpose(squeeze(curAvail));
    
    % Calculate per-axon probabilities
    curProbs = mnpdf(classConn, curAvail);
    
    % Calculate per-synapse probabilities
    curProbs = curProbs .^ (1 ./ synCounts);
    distProbs(:, curDistIdx) = curProbs;
end

%% Plotting
excProbs = distProbs(axonClasses(1).axonIds, :);
inhProbs = distProbs(axonClasses(2).axonIds, :);

fig = figure();
ax = axes(fig);
hold(ax, 'on');

plot(avail.dists, median(excProbs, 1), 'LineWidth', 2);
plot(avail.dists, median(inhProbs, 1), 'LineWidth', 2);

ylim(ax, [0, 1]);