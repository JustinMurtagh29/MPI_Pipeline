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

% Determine relative availabilities of target classes
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

% Count synapses per axon
% NOTE(amotta): Depending on the crazy stuff we're doing (e.g., excluding
% 'OtherDendrites') this value might be different from what's stored in
% conn.axonMeta.synCount.
synCounts = sum(classConn, 2);

%% Plot R² per axon and dendrite class
plotAxonClasses = 1:2;
plotTargetClasses = targetClasses;

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1100, 850];

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
    hold(curAx, 'on');
    
    for curTargetClassIdx = 1:numel(plotTargetClasses)
        plot( ...
            curAx, avail.dists / 1E3, ...
            curRsq(:, curTargetClassIdx), ...
            'LineWidth', 2);
    end
    
    ylim(curAx, [-0.2, 0.5]);
    xlabel('Radius (µm)');
    
    title( ...
        curAx, curAxonClass.title, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

ylabel(fig.Children(end), 'R²');
legend(fig.Children(1), ...
    arrayfun( ...
        @char, plotTargetClasses, ...
        'UniformOutput', false), ...
	'Location', 'NorthEast');

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'String', { ...
        info.filename; ...
        info.git_repos{1}.hash});