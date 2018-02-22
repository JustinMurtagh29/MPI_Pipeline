% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
availFile = '/tmpscratch/amotta/l4/2018-02-02-surface-availability-connectome-axons-18-a/axon-avail-data.mat';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a.mat');

minSynPre = 10;
info = Util.runInfo();

%% loading data
avail = load(availFile);
conn = load(connFile);

%% build class connectome
% try to replicate cass connectome
connectome = conn.connectome;

% count synapses
connectome.synCount = cellfun(@numel, connectome.synIdx);
connectome.synIdx = [];

% add target class to connectome
[~, denMetaRow] = ismember(connectome.edges(:, 2), conn.denMeta.id);
connectome.targetClass = conn.denMeta.targetClass(denMetaRow);

targetClasses = unique(conn.denMeta.targetClass);
[~, connectome.targetClassId] = ...
    ismember(connectome.targetClass, targetClasses);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);

%% build availabilities
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

%% likelihoods
axonTargetLik = cellfun(@numel, {avail.dists, conn.axons, targetClasses});
axonTargetLik = nan(axonTargetLik);

for curAxonId = 1:numel(conn.axons)
    curTargetSyns = classConnectome(curAxonId, :);
    curSynCount = conn.axonMeta.synCount(curAxonId);
    curLambdas = curSynCount * availabilities(:, :, curAxonId);
    
    for curClassIdx = 1:numel(targetClasses)
        axonTargetLik(:, curAxonId, curClassIdx) = poisspdf( ...
            curTargetSyns(curClassIdx), curLambdas(curClassIdx, :));
    end
end

%% build axon classes
conn.axonMeta.spineSynFrac = ...
    conn.axonMeta.spineSynCount ...
    ./ conn.axonMeta.synCount;

axonClasses = struct;
axonClasses(1).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac > 0.7);
axonClasses(1).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at least 70 %% onto spines'], minSynPre);

axonClasses(2).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac > 0.5);
axonClasses(2).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at least 50 %% onto spines'], minSynPre);

axonClasses(3).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac < 0.5);
axonClasses(3).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at most 50 %% onto spines'], minSynPre);
    
axonClasses(4).axonIds = find( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.spineSynFrac < 0.3);
axonClasses(4).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at most 30 %% onto spines'], minSynPre);

%% plot likelihoods
fig = figure();
for curTargetClassIdx = fliplr(1:numel(targetClasses))
    curTargetClass = targetClasses(curTargetClassIdx);
    
    curAx = subplot(1, numel(targetClasses), curTargetClassIdx);
    hold(curAx, 'on');
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonIds = axonClasses(curAxonClassIdx).axonIds;
        curMedianLik = axonTargetLik(:, curAxonIds, curTargetClassIdx);
        curMedianLik = median(curMedianLik, 2);
        
        plot(curAx, avail.dists / 1E3, curMedianLik, 'LineWidth', 2);
    end
    
    ylim(curAx, [0, 1]);
    title(curAx, char(curTargetClass));
    
    xlabel(curAx, 'r_{pred} (µm)');
    
    if curTargetClassIdx ~= 1
        yticklabels(curAx, {});
    end
    
    curAx.TickDir = 'out';
end

ylabel(curAx, { ...
    'Likelihood of target synapses';
    'under availability-based Poisson model'});
legend(curAx, {axonClasses.title}, 'Location', 'NorthWest');
annotation( ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% plot correlations
fig = figure();
for curTargetClassIdx = fliplr(1:numel(targetClasses))
    curTargetClass = targetClasses(curTargetClassIdx);
    
    curAx = subplot(1, numel(targetClasses), curTargetClassIdx);
    hold(curAx, 'on');
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonIds = axonClasses(curAxonClassIdx).axonIds;
        
        curAvails = availabilities(curTargetClassIdx, :, curAxonIds);
        curAvails = transpose(squeeze(curAvails));
        
        curSpecs = classConnectome(curAxonIds, :);
        curSpecs = curSpecs ./ sum(curSpecs, 2);
        
        curCorrVec = nan(numel(avail.dists), 1);
        for curDistIdx = 1:numel(avail.dists)
            curCorrMat = corrcoef( ...
                curAvails(:, curDistIdx), ...
                curSpecs(:, curTargetClassIdx));
            curCorrVec(curDistIdx) = curCorrMat(2);
        end
        
        plot(curAx, avail.dists / 1E3, curCorrVec, 'LineWidth', 2);
    end
    
    plot(curAx, avail.dists([1, end]) / 1E3, [0, 0], 'k--');
    
    title(curAx, char(curTargetClass));
    xlabel(curAx, 'r_{pred} (µm)');
    
    if curTargetClass == targetClasses(1)
        ylabel(curAx, { ...
            'Correlation between'; ...
            'availability and specificity'});
    else
        curAx.YTickLabel = {};
    end
    
    curAx.TickDir = 'out';
end

yMin = min(arrayfun(@(a) a.YLim(1), fig.Children));
yMax = max(arrayfun(@(a) a.YLim(2), fig.Children));
[fig.Children.YLim] = deal([yMin, yMax]);

legend(curAx, {axonClasses.title}, 'Location', 'NorthWest');
annotation( ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');