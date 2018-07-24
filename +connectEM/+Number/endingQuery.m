% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
skelDir = '/tmpscratch/amotta/l4/2018-07-23-ending-query-evaluation/nml-files';
nmlEvalFile = '/tmpscratch/amotta/l4/2018-07-23-ending-query-evaluation/all-ending-queries_v1.mat';

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
    voxelSize = param.raw.voxelSize;

    %% Find nml Files
    Util.log('Finding NML files');
    
    skelDirs = dir(skelDir);
    skelDirs = {skelDirs([skelDirs.isdir]).name};
    skelDirs(startsWith(skelDirs, '.')) = [];
    skelDirs = fullfile(skelDir, skelDirs);
    skelDirs = reshape(skelDirs, [], 1);

    nmlFiles = cellfun( ...
        @(skelDir) dir(fullfile(skelDir, '*.nml')), ...
        skelDirs, 'UniformOutput', false);
    nmlFiles = arrayfun( ...
        @(f) fullfile(f.folder, f.name), ...
        cat(1, nmlFiles{:}), 'UniformOutput', false);

    %% Processing
    Util.log('Processing NML files');

    nmlT = {'filePath', 'error', 'pathLenNm', 'timeMs'};
    nmlT = cell2struct(cell(numel(nmlFiles), numel(nmlT)), nmlT, 2);
   [nmlT.filePath] = deal(nmlFiles{:});

    parfor curIdx = 1:numel(nmlT)
        curFilePath = nmlT(curIdx).filePath;

        try
           [nmlT(curIdx).pathLenNm, nmlT(curIdx).timeMs] = ...
                forNmlFile(curFilePath, voxelSize);
        catch processingError
            nmlT(curIdx).error = processingError;
            nmlT(curIdx).pathLenNm = nan;
            nmlT(curIdx).timeMs = nan;
        end
    end
    
    nmlT = struct2table(nmlT, 'AsArray', true);

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
fprintf( ...
    'Time of flight: %g ± %g s (median = %g)\n', ...
    mean(nmlT.time), std(nmlT.time), median(nmlT.time));
fprintf( ...
    'Path length: %g ± %g µm (median = %g)\n', ...
    mean(nmlT.pathLen), std(nmlT.pathLen), median(nmlT.pathLen));
fprintf( ...
    'Flight speed: %g ± %g µm / s (median = %g)\n', ...
    mean(nmlT.speed), std(nmlT.speed), median(nmlT.speed));

% Per-project evaluation
[~, nmlT.projectName] = cellfun( ...
    @(path) fileparts(fileparts(path)), ...
    nmlT.filePath, 'UniformOutput', false);

projectT = table;
[projectT.name, ~, uniRows] = unique(nmlT.projectName);
projectT.taskCount = accumarray(uniRows, 1);
projectT.totalTimeH = accumarray(uniRows, nmlT.time) / 3600;
projectT.timePerTaskS = 3600 * projectT.totalTimeH ./ projectT.taskCount;

fprintf('\n');
disp(projectT)

fprintf('\n');
fprintf('Total number of tasks: %d\n', height(nmlT));
fprintf('Total number of work hours: %g h\n', sum(nmlT.time) / 3600);

%% Utility
function [lenNm, timeMs] = forNmlFile(path, voxelSize)
    nml = slurpNml(path);
    nodes = NML.buildNodeTable(nml);
    nodes.coord = nodes.coord .* voxelSize;
    
    % NOTE(amotta): Estimate time spent on this particular tracing. If
    % there are periods of inactivity that are 60 seconds or longer, they
    % are considered as separators between separate tracing sessions, and
    % temporarily pause the timer.
    %     Here, we're only considering flight mode tracings, where nodes
    % are automatically placed at a fixed (spatial) interval. A typical
    % inter-node interval is ~1 second. The above threshold should thus be
    % relatively robust.
    timeMs = diff(sort(nodes.time));
    timeMs(timeMs >= 60E3) = 0;
    timeMs = sum(timeMs);
    
    lenNm = cell2mat(cellfun( ...
        @(e) horzcat(e.source, e.target), ...
        nml.things.edges, 'UniformOutput', false));
   [~, lenNm] = ismember(lenNm, nodes.id);
   
    lenNm = ...
        nodes.coord(lenNm(:, 1), :) ...
      - nodes.coord(lenNm(:, 2), :);
    lenNm = sqrt(sum(lenNm .* lenNm, 2));
    
    % NOTE(amotta): For a period of time there was a bug in webKNOSSOS'
    % dynamic task switching, which occasionally introduced an edge when
    % the user was teleported to the location of the next task.
    %   To prevent this bug from influencing that path length measurement,
    % we're only counting edges below 10 µm. Since flight mode nodes are
    % automatically placed at a regular (sub-micron) spatial interval, this
    % heuristic should be fair.
    lenNm(lenNm >= 10E3) = 0;
    lenNm = sum(lenNm);
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
