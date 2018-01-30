% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';

connFile = fullfile( ...
    rootDir, 'connectomeState', ...
    'connectome_axons_18_a_with_den_meta.mat');

minSynPre = 10;
info = Util.runInfo();

%% loading data
conn = load(connFile);

%% preparing data
axonMask = (conn.axonMeta.synCount >= minSynPre);
classConnectome = conn.classConnectome(axonMask, :);

classNames = conn.denClasses;
classCount = numel(classNames);

classNames{strcmpi(classNames, 'whole cells')} = 'WCs';
classNames{strcmpi(classNames, 'apical dendrites')} = 'ADs';
classNames{strcmpi(classNames, 'smooth dendrite')} = 'SDs';

%% plotting
fig = figure;

for curRow = 1:classCount
    curRowSpecificities = ...
        classConnectome(:, curRow) ...
        ./ sum(classConnectome, 2);
    
    % virtually remove current class
    curClassConnectome = classConnectome;
    curClassConnectome(:, curRow) = 0;
    
    % calculate conditional specificity
    curCondSpecificities = ...
        curClassConnectome ...
        ./ sum(curClassConnectome, 2);
    
    for curCol = 1:classCount
        ax = subplot( ...
            classCount, classCount, ...
            curCol + (curRow - 1) * classCount);
        
        if curCol == curRow
            ax.Visible = 'off';
            continue;
        end
        
        scatter( ...
            ax, curRowSpecificities, ...
            curCondSpecificities(:, curCol), ...
            '.');
        
        if curCol == 1 && curRow == classCount
            xticks(ax, [0, 1]);
            yticks(ax, [0, 1]);
        end
        
        xlim(ax, [0, 1]);
        ylim(ax, [0, 1]);
        
        xlabel(sprintf('S(%s)', classNames{curRow}));
        ylabel(sprintf('S_{rem}(%s)', classNames{curCol}));
        
        if curCol == 1 && curRow == classCount
            xticks(ax, [0, 1]);
            yticks(ax, [0, 1]);
        else
            xticks(ax, []);
            yticks(ax, []);
        end
    end
end

annotation(...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {'Conditional target specificities'; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center')

fig.Position(3:4) = [1117, 1117];

%% plot proposed by Emmanuel
fig = figure;
ax = axes(fig);

for curIdx = 1:classCount
    curName = classNames{curIdx};
    
    % specificities for `x` class
    curSpecs = classConnectome(:, curIdx) ./ sum(classConnectome, 2);

    % specificities for remaining synapses
    curRemSpecs = classConnectome;
    curRemSpecs(:, curIdx) = [];
    
    % remove axons that have no other synapses
    % TODO(amotta): discuss if this is allowed
    curSpecs(~any(curRemSpecs, 2)) = [];
    curRemSpecs(~any(curRemSpecs, 2), :) = [];
    
    % calculate specificity
    curRemSpecs = curRemSpecs ./ sum(curRemSpecs, 2);

    curRemNames = classNames;
    curRemNames(curIdx) = [];

    % sort by increasing specificity for `x` class
    [~, sortIds] = sort(curSpecs, 'descend');
    curSpecs = curSpecs(sortIds);
    curRemSpecs = curRemSpecs(sortIds, :);

    ySpecsRed = cell2mat(arrayfun( ...
        @(i) median(curRemSpecs(1:i, :), 1), ...
        reshape(1:numel(curSpecs), [], 1), ...
        'UniformOutput', false));
    
    ax = subplot(1, classCount, curIdx);
    hold(ax, 'on');
    
    curColors = ax.ColorOrder;
    curColors(curIdx, :) = [];
    
    for curPlotIdx = 1:size(ySpecsRed, 2)
        plot( ...
            ax, ...
            curSpecs, ySpecsRed(:, curPlotIdx), ...
            'Color', curColors(curPlotIdx, :), ...
            'Marker', 'square', ...
            'MarkerSize', 3);
    end
    
    legend(curRemNames, 'Location', 'East');

    xlim(ax, [0, 1]); xticks(ax, [0, 1]);
    xlabel(sprintf('%s specificity threshold', curName));
    
    ylim(ax, [0, 1]); yticks(ax, [0, 1]);
    ax.TickDir = 'out';
    
    if curIdx == 1
        ylabel('Median specificity of remaining synapses');
    end
end

annotation(...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {
        'Thresholded conditional target specificities';
        info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

fig.Position(3:4) = [1916, 522];