% This script is based on 
%   +connectEM/+Connectome/plotGeometricPredictability.m
%   commit 8b9c96a8b57181c33f73cd411cf908138a76ca89
%
%   +connectEM/+Figure/coinnervation.m
%   commit f2418a0f552da11ce154faeffd01d65a4d867a11
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
% From +connectEM/calibrateNeuritePathLength.m (4397537e608a53268ab940e500551b3ad7cc33f2)
lengthFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
% From +connectEM/+Connectome/plotGeometricPredictability.m (1deca56c06f09ee8947e5cb9163dca09a124a4c1)
areaFile = '/tmpscratch/amotta/l4/2018-07-18-surface-availability-for-connectome-v7-partially-split/axon-availability_v2.mat';

[~, lengthFile] = fileparts(lengthFile);
lengthFile = sprintf('%s_pathLengths.mat', lengthFile);
lengthFile = fullfile(fileparts(connFile), lengthFile);

% NOTE(amotta): See
% connectEM.calibrateNeuritePathLength
% 107ae8f7e3819ae6920a8c7b19b4895868bbe476
dendCorrCoeff = 0.454;

minSynPre = 10;

synTypes = { ...
    'All', ...
    'PrimarySpine', ...
    'SecondarySpine', ...
    'Shaft', ...
    'Soma'};

axonClasses = { ...
    'Corticocortical', 'CC';
    'Thalamocortical', 'TC';
    'Inhibitory', 'Inh';
    'Other', 'Other'};

targetClasses = { ...
    'Somata', 'SO';
    'ProximalDendrite', 'PD'; ...
    'SmoothDendrite', 'SD'; ...
    'ApicalDendrite', 'AD'; ...
    'AxonInitialSegment', 'AIS'; ...
    'OtherDendrite', 'Other'};

axonTags = axonClasses(:, 2);
axonClasses = axonClasses(:, 1);

targetTags = targetClasses(:, 2);
targetClasses = targetClasses(:, 1);

% Utility
renorm = @(v) v / sum(v);

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = ...
    connectEM.Connectome.load(param, connFile);
conn = ...
    connectEM.Connectome.prepareForSpecificityAnalysis( ...
        conn, [], 'minSynPre', minSynPre);

lengths = load(lengthFile);
areas = load(areaFile);

%% Prepare availabilities
clear cur*;

[~, curIds] = ismember(targetClasses, areas.targetClasses);
postSynAreaFracs = areas.axonAvail(curIds, end, 1);
postSynAreaFracs = renorm(postSynAreaFracs(:));

% NOTE(amotta): The path lengths of soma-seeded agglomerates reported in
% the paper were derived from the ground truth skeleton tracings. Here,
% we're using the automatically calculated values (with calibration).
assert(isequal(numel(conn.dendrites), numel(lengths.dendritePathLengths)));

postSynLengths = lengths.dendritePathLengths;
[curMask, curIds] = ismember(conn.denMeta.targetClass, targetClasses);
postSynLengths = accumarray(curIds(curMask), postSynLengths(curMask));
postSynLengths(strcmpi(targetClasses, 'Somata')) = 0;
postSynLengthFracs = renorm(postSynLengths);

%% Run analysis
clear cur*;

for curSynIdx = 1:numel(synTypes)
    curSynType = synTypes{curSynIdx};
    curSynLut = syn.synapses.type == curSynType;
    if strcmpi(curSynType, 'all'); curSynLut(:) = true; end
    
    curConn = conn;
    curConn.connectome.synIdx = cellfun( ...
        @(ids) ids(curSynLut(ids)), ...
        curConn.connectome.synIdx, ...
        'UniformOutput', false);
    
    curClassConn = ...
        connectEM.Connectome.buildClassConnectome( ...
            curConn, ...
            'axonClasses', axonClasses, ...
            'targetClasses', targetClasses);
        
    curSynCount = sum(curClassConn(:));
    curTitleStem = sprintf('%s (n = %d)', curSynType, curSynCount);
        
    %% Synapse frequencies
    curTitle = strcat(curTitleStem, ' versus synapse type frequencies');
    
    curAxonFracs = renorm(sum(curClassConn, 2));
    curTargetFracs = renorm(sum(curClassConn, 1));

    curRelClassConn = curClassConn ./ sum(curClassConn(:));
    curExpClassConn = curAxonFracs .* curTargetFracs;
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
        
    plotIt( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
    
    %% Membrane frequencies
    curTitle = strcat(curTitleStem, ' versus membrane contributions');
    
    curAxonFracs = renorm(sum(curClassConn, 2));
    curTargetFracs = reshape(postSynAreaFracs, 1, []);
    
    curRelClassConn = curClassConn ./ sum(curClassConn, 2);
    curExpClassConn = repmat(curTargetFracs, size(curClassConn, 1), 1);
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
    
    plotIt( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
    
    %% Path length frequencies
    curTitle = strcat(curTitleStem, ' versus path length contributions');
    
    curAxonFracs = renorm(sum(curClassConn, 2));
    curTargetFracs = reshape(postSynLengthFracs, 1, []);
    
    curRelClassConn = curClassConn ./ sum(curClassConn, 2);
    curExpClassConn = repmat(curTargetFracs, size(curClassConn, 1), 1);
    
    curCorrCoeffs = curRelClassConn ./ curExpClassConn;
    curCorrCoeffs(curRelClassConn == 0 & curExpClassConn == 0) = 1;
    
    plotIt( ...
        info, curTitle, ...
        axonTags, curAxonFracs, ...
        targetTags, curTargetFracs, ...
        curRelClassConn, curExpClassConn, curCorrCoeffs);
end

%% Plotting
function plotIt( ...
        info, titleStr, ...
        axonTags, axonFracs, ...
        targetTags, targetFracs, ...
        relClassConn, expClassConn, corrCoeffs)
    rows = numel(axonTags);
    cols = numel(targetTags);
    frac = rows / cols;

    colorMat = log10(corrCoeffs);
    colorMap = connectEM.Figure.redBlue(129);
    colorLim = 1;
    % colorLim = ceil(max(abs(colorMat(isfinite(colorMat)))));

    fig = figure();
    ax = axes(fig);

    imshow( ...
        colorMat, [-colorLim, +colorLim], ...
        'Colormap', colorMap, ...
        'Parent', ax);

    fig.Color = 'white';
    fig.Position(3:4) = 1000 .* [1, frac];

    ax.Visible = 'on';
    ax.TickDir = 'out';
    ax.Box = 'off';

    ax.XAxisLocation = 'top';
    ax.XTick = 1:size(colorMat, 2);
    ax.XTickLabel = arrayfun( ...
        @(t, p) sprintf('%s (%.2g)', t{1}, 100 * p), ...
        targetTags(:), targetFracs(:), 'UniformOutput', false);
    ax.XTickLabelRotation = 90;

    ax.YTick = 1:size(colorMat, 1);
    ax.YTickLabel = arrayfun( ...
        @(t, p) sprintf('%s (%.2g)', t{1}, 100 * p), ...
        axonTags(:), axonFracs(:), 'UniformOutput', false);
    ax.Position = [0.1, 0.01, 0.75, 0.75];

    for curIdx = 1:numel(colorMat)
       [curRow, curCol] = ind2sub(size(colorMat), curIdx);
        curEdgeColor = 'none';

        curBoxSize = ax.Position(3:4) ./ [cols, rows];
        curOff = [curCol, numel(axonTags) - curRow + 1];
        curOff = ax.Position(1:2) + (curOff - 1) .* curBoxSize;

        curAnn = annotation(fig, ...
            'textbox', [curOff, curBoxSize], ...
            'String', { ...
                sprintf('%.2g', 100 * relClassConn(curIdx)); ...
                sprintf('(%.2g)', 100 * expClassConn(curIdx))});

        curAnn.HorizontalAlignment = 'center';
        curAnn.VerticalAlignment = 'middle';
        curAnn.EdgeColor = curEdgeColor;
        curAnn.Color = 'black';
        curAnn.FontSize = 12;
        curAnn.LineWidth = 2;
    end

    cbar = colorbar('peer', ax);
    cbar.Label.String = '"Affinity" (measured / expected)';
    cbar.Ticks = (-colorLim):colorLim;
    cbar.TickLabels = arrayfun( ...
        @(f) sprintf('%g', 10 ^ f), ...
        cbar.Ticks, 'UniformOutput', false);
    cbar.TickDirection = 'out';
    cbar.Position = [0.86, 0.1, 0.02, 0.8];
    cbar.Position([2, 4]) = ax.Position([2, 4]);

    title(ax, ...
        {info.filename; info.git_repos{1}.hash; titleStr}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end
