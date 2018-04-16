% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons_18_a_ax_spine_syn_clust.mat');
synFile = fullfile(rootDir, 'connectomeState', 'SynapseAgglos_v3_ax_spine_clustered_classified.mat');

targetClasses = { ...
    'Somata', 'SO'; ...
    'WholeCell', 'WC'; ...
    'ApicalDendrite', 'AD'; ...
    'SmoothDendrite', 'SD'; ...
    'AxonInitialSegment', 'AIS'};
targetLabels = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile, synFile);

%% Plotting
for curIdx = 1:numel(axonClasses)
    curAxonClass = axonClasses(curIdx);
    curAxonIds = curAxonClass.axonIds;
    curSpecs = curAxonClass.specs;
    
    curIds = fieldnames(curSpecs);
   [~, curIds] = ismember(curIds, targetClasses);
   
    curIds = nonzeros(curIds);
    curSpecClasses = targetClasses(curIds);
    curSpecLabels = targetLabels(curIds);
    
    curSpecAxonIds = cellfun( ...
        @(name) curSpecs.(name).axonIds, ...
        curSpecClasses, 'UniformOutput', false);
    
    curSpecMat = 1:numel(curSpecLabels);
   [curA, curB] = ndgrid(curSpecMat, curSpecMat);
    
    curSpecMat = zeros(1 + numel(curSpecLabels));
    curSpecMat(1:(end - 1), 1:(end - 1)) = cellfun( ...
        @(idsOne, idsTwo) numel(intersect(idsOne, idsTwo)), ...
        curSpecAxonIds(curA), curSpecAxonIds(curB));
    
    curNonSpecAxonIds = setdiff(curAxonIds, cell2mat(curSpecAxonIds));
    curSpecMat(end, end) = numel(curNonSpecAxonIds);
    
    plotIt(info, curAxonClass, curSpecLabels, curSpecMat)
end

%% Utility
function plotIt(info, axonClass, specLabels, specMat)
    axonClassTitle = axonClass.title;
    axonCount = numel(axonClass.axonIds);
    specLabels{end + 1} = 'None';
    
    %% Plot matrix
    fig = figure();
    ax = axes(fig);
    
    specMatRel = specMat / axonCount;
    imshow(specMatRel, 'Parent', ax, 'Colormap', parula);

    fig.Color = 'white';
    fig.Position(3:4) = [660, 660];

    ax.Visible = 'on';
    ax.Box = 'off';
    ax.TickDir = 'out';
    
    ax.XAxisLocation = 'top';
    ax.XTick = 1:numel(specLabels);
    ax.XTickLabel = specLabels;
    ax.XTickLabelRotation = 90;

    ax.YTick = 1:numel(specLabels);
    ax.YTickLabel = specLabels;

    axis(ax, 'square');
    ax.Position = [0.1, 0.05, 0.8, 0.8];
    
    for curIdx = 1:numel(specMat)
       [curRow, curCol] = ind2sub(size(specMat), curIdx);

        curBoxSize = ax.Position(3:4) / numel(specLabels);
        curOff = [curCol, numel(specLabels) - curRow + 1];
        curOff = ax.Position(1:2) + (curOff - 1) .* curBoxSize;

        annotation(fig, ...
            'textbox', [curOff, curBoxSize], ...
            'String', { ...
                sprintf('%d', specMat(curIdx)); ...
                sprintf('%.1f %%', 100 * specMatRel(curIdx))}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'EdgeColor', 'none', ...
            'Color', 'white', ...
            'FontSize', 12);
    end
    
    cbar = colorbar('peer', ax);
    cbar.TickDirection = 'out';
    cbar.Position = [0.925, 0.05, 0.02, 0.8];
    
    title(ax, ...
        {info.filename; info.git_repos{1}.hash; axonClassTitle}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end