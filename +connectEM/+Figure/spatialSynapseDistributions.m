% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

binEdges = linspace(-50, 50, 51);

synTypes = { ...
    'Corticocortical', ...
    'Thalamocortical', ...
    'Inhibitory'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

%% Preparing data
synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.synType = conn.axonMeta.axonClass(synT.preAggloId);

[~, synT.synType] = ismember(synT.synType, synTypes);
synT(~synT.synType, :) = [];

synPos = connectEM.Synapse.calculatePositions(param, syn);
synT.pos = synPos(synT.id, :);
clear synPos;

% Physical units, relative to center
synT.pos = synT.pos - mean(param.bbox, 2)';
synT.pos = synT.pos .* param.raw.voxelSize ./ 1E3;

%% Plotting data
dimLabels = {'X', 'Y', 'Z'};

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1500, 675];

plotHist = @(ax, data) ...
    histogram(ax, ...
        'BinEdges', binEdges, ...
        'BinCounts', data, ...
        'Orientation', 'horizontal', ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    
for curDim = 1:3
    curPosId = synT.pos(:, curDim);
    curPosId = discretize(curPosId, binEdges);
    
    curSynCounts = accumarray( ...
        cat(2, curPosId, synT.synType), ...
        1, [numel(binEdges) - 1, numel(synTypes)]);
    curSynRatios = [ ...
        curSynCounts(:, 3) ./ sum(curSynCounts, 2), ...
        curSynCounts(:, 2) ./ sum(curSynCounts(:, 1:2), 2)];
    
    % Synapse probabilities
    curAx = subplot(2, 3, curDim);
    curAx.YAxis.Direction = 'reverse';
    hold(curAx, 'on');
    
    for curTypeId = 1:size(curSynCounts, 2)
        curProbs = curSynCounts(:, curTypeId);
        curProbs = curProbs / sum(curProbs);
        plotHist(curAx, curProbs);
    end
    
    % Synapse ratios
    curAx = subplot(2, 3, curDim + 3);
    curAx.YAxis.Direction = 'reverse';
    hold(curAx, 'on');
    
    for curRatioId = 1:size(curSynRatios, 2)
        curRatios = curSynRatios(:, curRatioId);
        curRatios(isnan(curRatios)) = 0;
        plotHist(curAx, curRatios);
    end
    
    ylabel(curAx, sprintf('%s (µm)', dimLabels{curDim}));
end

axes = flip(fig.Children);

% Make square
[axes.PlotBoxAspectRatio] = deal([1, 1, 1]);
[axes.DataAspectRatioMode] = deal('auto');

% X labels
xlabel(axes(1), 'Probability');
xlabel(axes(2), 'Ratio');

[axes(1:2:end).XLim] = deal([0, 0.06]);
[axes(2:2:end).XLim] = deal([0, 0.4]);

% Legends
curLeg = addLegend(axes(end - 1), synTypes, 'Location', 'EastOutside');
curLeg.Box = 'off';

curLeg = addLegend(axes(end), ...
   {'Inh / (Inh + Exc)', 'TC / (TC + CC)'}, ...
    'Location', 'EastOutside');
curLeg.Box = 'off';

annotation( ...
    fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
	'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Utilities
function leg = addLegend(ax, legends, varargin)
    axPos = ax.Position;
    leg = legend(ax, legends, varargin{:});
    ax.Position = axPos;
end
