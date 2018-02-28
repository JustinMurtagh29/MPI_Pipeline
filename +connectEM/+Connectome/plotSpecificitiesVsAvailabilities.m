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

% renumber classes alphabetically
conn.denMeta.targetClass = reordercats(conn.denMeta.targetClass);
targetClasses = unique(conn.denMeta.targetClass);

[~, connectome.targetClassId] = ...
    ismember(connectome.targetClass, targetClasses);

classConnectome = accumarray( ...
    cat(2, connectome.edges(:, 1), connectome.targetClassId), ...
    connectome.synCount, [numel(conn.axons), numel(targetClasses)]);
specificities = classConnectome ./ sum(classConnectome, 2);

%% build availabilities
[~, classIds] = ismember(targetClasses, avail.targetClasses);
availabilities = avail.axonAvail(classIds, :, :);
availabilities = availabilities ./ sum(availabilities, 1);

%% plot availability vs. distance
axonIds = find(conn.axonMeta.synCount >= 10);
cmap = parula(101);

fig = figure();
fig.Color = 'white';

for classIdx = 1:numel(targetClasses)
    className = char(targetClasses(classIdx));
    
   [~, sortIds] = sort(specificities(axonIds, classIdx), 'ascend');
    axonIds = axonIds(sortIds);
    
    ax = subplot(1, numel(targetClasses), classIdx);
    
    ax.Color = 'black';
    ax.TickDir = 'out';
    
    colormap(ax, cmap);
    hold(ax, 'on');

    for curAxonId = reshape(axonIds, 1, [])
        curAvails = shiftdim(availabilities(classIdx, :, curAxonId));

        curColor = specificities(curAxonId, classIdx);
        curColor = cmap(round(100 * curColor) + 1, :);

        plot(ax, avail.dists / 1E3, curAvails, 'Color', curColor);
    end

    xlabel(ax, {className; 'r_{pred} (µm)'});
    xlim(ax, [0, avail.dists(end)] / 1E3);
end

% Y axis
yMax = max(arrayfun(@(a) a.YLim(end), fig.Children));
[fig.Children.YLim] = deal([0, yMax]);
ylabel(fig.Children(end), 'Availability');

% colorbar
oldPos = ax.Position;
cbar = colorbar(ax);
cbar.Label.String = 'Specificity';
ax.Position = oldPos;

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% scatter plot with specificity vs. availability
axonIds = find(conn.axonMeta.synCount >= 10);
plotDists = [1E3 .* [1, 10, 20, 50], avail.dists(end)];
[~, plotDists] = ismember(plotDists, avail.dists);

classCount = numel(targetClasses);
distCount = numel(plotDists);

fig = figure();
fig.Color = 'white';

for curDistIdx = 1:distCount
    curDistId = plotDists(curDistIdx);
    curDistUm = avail.dists(curDistId) / 1E3;
    
    for curClassIdx = 1:classCount
        curAvail = squeeze(availabilities( ...
            curClassIdx, curDistId, axonIds));
        curSpecs = specificities(axonIds, curClassIdx);
        
        curFit = fit( ...
            curSpecs, curAvail, 'poly1', ...
            'Weights', conn.axonMeta.synCount(axonIds));
        
        curAxIdx = curClassIdx + (curDistIdx - 1) * classCount;
        curAx = subplot(distCount, classCount, curAxIdx);
        curAx.TickDir = 'out';
        
        hold(curAx, 'on');
        xlim(curAx, [0, 1]);
        ylim(curAx, [0, 1]);
        
        scatter(curAx, curSpecs(:), curAvail(:), '.');
        plot(curAx, xlim(curAx), curFit(xlim(curAx)));
        annotation( ...
            fig, ...
            'textbox', curAx.Position, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
            'String', sprintf('y = %.2fx + %.2f', curFit.p1, curFit.p2));
        
        if curClassIdx == 1
            ylabel(curAx, {'Availability'; ...
                sprintf('r_{pred} = %d µm', curDistUm)});
        else
            yticklabels(curAx, {});
        end
        
        if curDistIdx == distCount
            xlabel(curAx, {'Specificity'; ...
                char(targetClasses(curClassIdx))});
        end
    end
end

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% plot fraction of specificity explained by availability
axonIds = find(conn.axonMeta.synCount >= 10);
plotDists = [1E3 .* [1, 10, 20, 50], avail.dists(end)];
[~, plotDists] = ismember(plotDists, avail.dists);

classCount = numel(targetClasses);
distCount = numel(plotDists);

maxX = nan(classCount, 1);
for curClassIdx = 1:classCount
    curAvail = squeeze(availabilities( ...
        curClassIdx, plotDists, axonIds))';
    curSpecs = specificities(axonIds, curClassIdx);
    
    curFracExp = reshape(curAvail ./ curSpecs, [], 1);
    curFracExp = curFracExp(~isinf(curFracExp));
    
    maxX(curClassIdx) = prctile(curFracExp, 99);
end

maxX = max(maxX, 1);

fig = figure();
fig.Color = 'white';

for curDistIdx = 1:distCount
    curDistId = plotDists(curDistIdx);
    curDistUm = avail.dists(curDistId) / 1E3;
    
    for curClassIdx = 1:classCount
        curAvail = squeeze(availabilities( ...
            curClassIdx, curDistId, axonIds));
        curSpecs = specificities(axonIds, curClassIdx);
        
        curFracExp = curAvail(:) ./ curSpecs(:);
        curMaxX = maxX(curClassIdx);
        
        curAxIdx = curClassIdx + (curDistIdx - 1) * classCount;
        curAx = subplot(distCount, classCount, curAxIdx);
        curAx.TickDir = 'out';
        
        hold(curAx, 'on');
        xlim(curAx, [0, curMaxX]);
        
        histogram( ...
            curAx, curFracExp, linspace(0, curMaxX, 21), ...
            'DisplayStyle', 'stairs', 'LineWidth', 2);
        plot(curAx, [1, 1], ylim(curAx), 'k--');
        
        if curClassIdx == 1
            ylabel(curAx, {sprintf( ...
                'r_{pred} = %d µm', curDistUm); 'Axons'});
        end
        
        switch curDistIdx
            case 1
                title(curAx, char(targetClasses(curClassIdx)));
            case distCount
                xlabel(curAx, 'Availability / Specificity');
        end
    end
end

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

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
  & conn.axonMeta.spineSynFrac < 0.3);
axonClasses(2).title = sprintf( ...
   ['axons with ≥ %d synapses and ', ...
    'at most 30 %% onto spines'], minSynPre);

%% plot cosine similarity
fig = figure();
for curTargetClassIdx = fliplr(1:numel(targetClasses))
    curTargetClass = targetClasses(curTargetClassIdx);
    
    curAx = subplot(1, numel(targetClasses), curTargetClassIdx);
    hold(curAx, 'on');
    
    for curAxonClassIdx = 1:numel(axonClasses)
        curAxonIds = axonClasses(curAxonClassIdx).axonIds;
        
        curAvail = squeeze(availabilities( ...
            curTargetClassIdx, :, curAxonIds));
        curSpecs = transpose(specificities( ...
            curAxonIds, curTargetClassIdx));
        
        curCos = 1 - pdist2(curSpecs, curAvail, 'cosine');
        plot(curAx, avail.dists / 1E3, curCos, 'LineWidth', 2);
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
    'Cosine similarity between';
    'availability and specificity'});
legend(curAx, {axonClasses.title}, 'Location', 'NorthWest');

annotation( ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');
