% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outFile = '/tmpscratch/amotta/l4/2018-05-08-ending-query-evaluation/axon-state_60.mat';

info = Util.runInfo();

%% Running analysis
if ~exist(outFile, 'file')
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

    Util.saveStruct(outFile, out);
    Util.protect(outFile);
else
    out = load(outFile);
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

% Plot time
ax = subplot(1, 3, 1);
axis(ax, 'square');
hold(ax, 'on');

binEdges = linspace(0, 100, 101);
meanTime = mean(nmlT.time);
medianTime = median(nmlT.time);

histogram( ...
    ax, nmlT.time, ...
    'BinEdges', binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
plot( ...
    ax, repelem(meanTime, 2, 1), ax.YLim, ...
    'Color', 'red', 'LineStyle', '--');
plot( ...
    ax, repelem(medianTime, 2, 1), ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Time of flight (s)');

yticklabels(arrayfun( ...
    @num2str, yticks(ax), ...
    'UniformOutput', false));
ylabel(ax, 'Queries');
ax.TickDir = 'out';

% Plot path length
ax = subplot(1, 3, 2);
axis(ax, 'square');
hold(ax, 'on');

binEdges = linspace(0, 20, 101);
meanPathLen = mean(nmlT.pathLen);
medianPathLen = median(nmlT.pathLen);

histogram( ...
    ax, nmlT.pathLen, ...
    'BinEdges', binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
plot( ...
    ax, repelem(meanPathLen, 2, 1), ax.YLim, ...
    'Color', 'red', 'LineStyle', '--');
plot( ...
    ax, repelem(medianPathLen, 2, 1), ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Length of flight (µm)');

yticklabels(arrayfun( ...
    @num2str, yticks(ax), ...
    'UniformOutput', false));
ylabel(ax, 'Queries');
ax.TickDir = 'out';

% Plot speed
ax = subplot(1, 3, 3);
axis(ax, 'square');
hold(ax, 'on');

binEdges = linspace(0, 2, 101);
meanSpeed = mean(nmlT.speed);
medianSpeed = median(nmlT.speed);

histogram( ...
    ax, nmlT.speed, ...
    'BinEdges', binEdges, ...
    'DisplayStyle', 'stairs', ...
    'LineWidth', 2);
hMean = plot( ...
    ax, repelem(meanSpeed, 2, 1), ax.YLim, ...
    'Color', 'red', 'LineStyle', '--');
hMedian = plot( ...
    ax, repelem(medianSpeed, 2, 1), ax.YLim, ...
    'Color', 'black', 'LineStyle', '--');

xlim(ax, binEdges([1, end]));
xlabel(ax, 'Flight speed (mm / h)');

yticklabels(arrayfun( ...
    @num2str, yticks(ax), ...
    'UniformOutput', false));
ylabel(ax, 'Queries');
ax.TickDir = 'out';

leg = legend( ...
    [hMean, hMedian], ...
    'Mean', 'Median', ...
    'Location', 'NorthEast');
leg.Box = 'off';

annotation(fig, ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'String', {info.filename; info.git_repos{1}.hash}, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center');

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

%% Plotting results
