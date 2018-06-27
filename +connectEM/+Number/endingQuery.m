% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
nmlEvalFile = '/tmpscratch/amotta/l4/2018-05-08-ending-query-evaluation/axon-state_60.mat';

% See https://gitlab.mpcdf.mpg.de/connectomics/pipeline/blob/c3e5dacf542e71c452b9b8d7e3fe5f63bd5b8e0c/+connectEM/setQueryState.m
axonsBeforeQueryFile = fullfile(rootDir, 'aggloState', 'axons_04.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Number of axon fragments at beginning
axonsBeforeQuery = load(axonsBeforeQueryFile, 'indBigAxons');
numFiveUmAxonFragments = sum(axonsBeforeQuery.indBigAxons) %#ok

%% Running analysis
if ~exist(nmlEvalFile, 'file')
    %% Loading data
    param = load(fullfile(rootDir, 'allParameter.mat'));
    param = param.p;

    %% Find nml Files
    Util.log('Finding NML files');

    skelDirs = connectEM.setQueryState('6.0');
    skelDirs = fullfile('/gaba', reshape(skelDirs, [], 1));

    nmlFiles = cellfun( ...
        @(skelDir) dir(fullfile(skelDir, '*.nml')), ...
        skelDirs, 'UniformOutput', false);
    nmlFiles = arrayfun( ...
        @(f) fullfile(f.folder, f.name), ...
        cat(1, nmlFiles{:}), 'UniformOutput', false);

    %% Processing
    Util.log('Processing NML files');

    nmlT = table;
    nmlT.filePath = nmlFiles;
    nmlT.pathLenNm = nan(height(nmlT), 1);
    nmlT.timeMs = nan(height(nmlT), 1);

    tic;
    for curIdx = 1:height(nmlT)
        curFilePath = nmlT.filePath{curIdx};

        try
           [nmlT.pathLenNm(curIdx), nmlT.timeMs(curIdx)] = ...
                forNmlFile(curFilePath, param.raw.voxelSize);
        catch
            warning('Error in %s', curFilePath);
        end

        Util.progressBar(curIdx, height(nmlT));
    end

    %% Save result
    out = struct;
    out.nmlT = nmlT;
    out.info = info;

    Util.saveStruct(nmlEvalFile, out);
    Util.protect(nmlEvalFile);
else
    out = load(nmlEvalFile);
end

%% Plotting results
nmlT = out.nmlT;

% Get rid of invalid tracings
nmlT(isnan(nmlT.pathLenNm) | nmlT.pathLenNm == 0, :) = [];
nmlT(isnan(nmlT.timeMs) | nmlT.timeMs == 0, :) = [];

nmlT.time = nmlT.timeMs / 1E3;
nmlT.pathLen = nmlT.pathLenNm / 1E3;
nmlT.speed = (nmlT.pathLen / 1E3) ./ (nmlT.time / 3600);

fig = figure();
fig.Color = 'white';
fig.Position(3:4) = [1200, 440];

% Plot time
ax = subplot(1, 3, 1);
binEdges = linspace(0, 100, 101);
plotData(ax, 'Time of flight', 's', nmlT.time, binEdges);

% Plot path length
ax = subplot(1, 3, 2);
binEdges = linspace(0, 15, 101);
plotData(ax, 'Length of flight', 'µm', nmlT.pathLen, binEdges);

% Plot speed
ax = subplot(1, 3, 3);
binEdges = linspace(0, 2, 101);
plotData(ax, 'Flight speed', 'mm / h', nmlT.speed, binEdges);

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', { ...
        info.filename; info.git_repos{1}.hash; ...
        sprintf('Evaluation of %d ending queries', height(out.nmlT))}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

%% Quantitative evaluation
meanQueryPathLengthUm = mean(nmlT.pathLen) %#ok
meanQuerySpeedUmPerS = mean(nmlT.pathLen ./ nmlT.time) %#ok
meanQueryTime = meanQueryPathLengthUm / meanQuerySpeedUmPerS %#ok

%% Utility
function [lenNm, timeMs] = forNmlFile(path, voxelSize)
    nml = slurpNml(path);
    nodes = NML.buildNodeTable(nml);
    nodes.coord = nodes.coord .* voxelSize;
    
    timeMs = sort(nodes.time);
    timeMs = timeMs(end) - timeMs(2);
    
    lenNm = cell2mat(cellfun( ...
        @(e) horzcat(e.source, e.target), ...
        nml.things.edges, 'UniformOutput', false));
   [~, lenNm] = ismember(lenNm, nodes.id);
   
    lenNm = ...
        nodes.coord(lenNm(:, 1)) ...
      - nodes.coord(lenNm(:, 2));
    lenNm = sum(lenNm .* lenNm, 2);
    lenNm = sum(sqrt(lenNm));
end

function plotData(ax, name, unit, values, binEdges)
    meanVal = mean(values);
    medianVal = median(values);
    
    axis(ax, 'square');
    hold(ax, 'on');

    histogram( ...
        ax, values, ...
        'BinEdges', binEdges, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1);
    meanPlot = plot( ...
        ax, repelem(meanVal, 2, 1), ax.YLim, ...
        'Color', 'red', 'LineStyle', '--');
    medianPlot = plot( ...
        ax, repelem(medianVal, 2, 1), ax.YLim, ...
        'Color', 'black', 'LineStyle', '--');

    xlim(ax, binEdges([1, end]));
    xlabel(ax, sprintf('%s (%s)', name, unit));
    
    yticklabels(arrayfun( ...
        @num2str, yticks(ax), ...
        'UniformOutput', false));
    ylabel(ax, 'Queries');
    ax.TickDir = 'out';
    
    leg = legend( ...
       [meanPlot, medianPlot], ...
        sprintf('Mean (%.2f %s)', meanVal, unit), ...
        sprintf('Median (%.2f %s)', medianVal, unit), ...
        'Location', 'NorthOutside');
    leg.Box = 'off';
end
