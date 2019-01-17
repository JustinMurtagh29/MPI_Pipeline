% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputMapFile = '/tmpscratch/amotta/l4/2018-07-26-tracing-based-output-maps/20180726T190355_results.mat';

binSizeUm = 20;

% Set output directory to save PNG and EPS files of plots
outputDir = '';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

outputMap = load(outputMapFile);
connFile = outputMap.info.param.connFile;
axonData = outputMap.axonData;

[conn, syn] = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

%% Ad-hoc patching axon-specific data
% * Add NML file ID (as proxy of cell ID) to synapse tables
% * Add target class (we will later use the WholeCell and ApicalDendrite
%   categories as proxies to find likely synapses onto L4 vs. L5 cells).
curDendClasses = [ ...
    categorical({'OtherDendrite'}); ...
    conn.denMeta.targetClass];

for curId = 1:numel(axonData)
    axonData(curId).synapses.nmlId(:) = curId;
    
    curTargetClass = axonData(curId).synapses.id;
   [~, curTargetClass] = ismember(curTargetClass, synT.id);
    
    % NOTE(amotta): Synapses were picked up regardless of the postsynaptic
    % process. A subset of these synapses are onto targets that are missing
    % from the connectome. This subset is not in the synapse table.
    %   For now, these synapses are treated as being onto OtherDendrite.
    curTargetClass(curTargetClass > 0) = ...
        synT.postAggloId(curTargetClass(curTargetClass > 0));
    curTargetClass = curDendClasses(1 + curTargetClass);
    
    axonData(curId).synapses.targetClass = curTargetClass;
end

%% Ignore interneuron axons for PLASS analysis
% NOTE(amotta): It turns out that of the two interneurons in the dataset
% (whole cells 17 and 22), only one has part of its axon in the EM volume.
% But said axon is mostly myelinated and doesn't make a single synapse.
% Still, let's remove these tracings for completeness.
clear cur*;

interNeuronGtNmlFiles = { ...
    '5a796fcd67000090172d94f1', ... % Whole cell 17. Ground truth tracing in +connectEM/+WholeCell/+Data/border-cells_axon-dendrites-split/5a796fcd67000090172d94f1.nml
    '5a796fcd67000060172d94e0'};    % Whole cell 22. Ground truth tracing in +connectEM/+WholeCell/+Data/border-cells_axon-dendrites-split/5a796fcd67000060172d94e0.nml

[~, curNmlFiles] = cellfun( ...
    @fileparts, {axonData.nmlFile}, 'UniformOutput', false);
axonData(ismember(curNmlFiles, interNeuronGtNmlFiles)) = [];
clear cur*;

%% Build grand-average axon
grandAvgAxon = struct;
grandAvgAxon.nmlFile = 'grand average';
grandAvgAxon.pathLength = nan;

grandAvgAxon.synapses = cat(1, axonData.synapses);
grandAvgAxon.synapses = sortrows(grandAvgAxon.synapses, 'somaDist');

%% Perform statistical test
clear cur*;
curSomaDist = grandAvgAxon.synapses.somaDist;
curIsSpineSyn = syn.isSpineSyn(grandAvgAxon.synapses.id);

[~, curPVal] = ttest2( ...
    curSomaDist(curIsSpineSyn), ...
    curSomaDist(~curIsSpineSyn)) %#ok

%% Determine bin edges
binEdges = max(grandAvgAxon.synapses.somaDist / 1E3);
binEdges = binSizeUm * ceil(binEdges / binSizeUm);
binEdges = 0:binSizeUm:binEdges;

axonCounts = arrayfun(...
    @(e) numel(unique( ...
        grandAvgAxon.synapses.nmlId( ...
            grandAvgAxon.synapses.somaDist >= 1E3 * e))), ...
	binEdges);
axonCounts = axonCounts(1:(end - 1));

%% Plot axon
clear cur*;
curData = cat(1, axonData(:), grandAvgAxon);
curData = reshape(curData, 1, []);

curNumDigits = ceil(log10(1 + numel(curData)));

for curIdx = 1:numel(curData)
    curAxonData = curData(curIdx);
    curSynapses = curAxonData.synapses;
    if isempty(curSynapses); continue; end
    
    curSynapses.somaDist = curSynapses.somaDist / 1E3;
    curSynapses.ontoSpine = syn.isSpineSyn(curSynapses.id);
   [~, curAxonName] = fileparts(curAxonData.nmlFile);
   
    curPlotData = @(varargin) plotData( ...
        binEdges, curSynapses.somaDist, varargin{:});
    curReport = curIdx == numel(curData);

    curFig = figure();
    curFig.Color = 'white';
    curFig.Position(3:4) = [600, 590];
    
    % Panel #1 - Spine vs. shaft synapses
    curAx = subplot(3, 1, 1);
    curPlotData( ...
        curAx, 1 + not(curSynapses.ontoSpine), ...
        {'Spine synapse', 'Non-spine synapse'}, ...
        'classColors', {'magenta', 'black'}, ...
        'report', curReport);
    
    xlabel(curAx, 'Axonal path length to soma (µm)');
    curAx.YLim(1) = 0;

    title(curAx, { ...
        info.filename; info.git_repos{1}.hash; ...
        sprintf('Axon %s', curAxonName)}, ...
        'FontWeight', 'normal', 'FontSize', 10);
    
    % Panel #2 - L4 vs. L5 synapses. The WholeCell and ApicalDendrite
    % target classes serve as proxies for likely L4 and L5 cells, resp.
    curAx = subplot(3, 1, 2);

   [~, curSynClasses] = ismember( ...
       curSynapses.targetClass, ...
       {'ProximalDendrite', 'ApicalDendrite'});
    curPlotData( ...
        curAx, curSynClasses, ...
        {'Onto PD', 'Onto AD'}, ...
        'report', curReport);
    curAx.YLim(1) = 0;
    
    % Panel #3 - Number of cells contributing to bin
    curAx = subplot(3, 1, 3);
    plotHist(curAx, binEdges, axonCounts, 'EdgeColor', 'black');
    
    curAx.YDir = 'reverse';
    curAx.XAxisLocation = 'top';
    
    curAx.XLim = binEdges([1, end]);
    curAx.XTickLabel = {''};
    curAx.TickDir = 'out';
    curAx.Box = 'off';
    
    % Make sure all axes have same length in plot
    curAxes = findobj(curFig.Children, 'Type', 'Axes');
    curAxesPos = cell2mat(get(curAxes, {'Position'}));
    curAxesPos(:, 3) = min(curAxesPos(:, 3));
    set(curAxes, {'Position'}, num2cell(curAxesPos, 2));
    set(curAxes, 'XLim', binEdges([1, end]));
    
    if ~isempty(outputDir)
        curFigFileName = strrep(curAxonName, ' ', '-');
        curFigFileName = fullfile(outputDir, sprintf( ...
            '%0*d_%s', curNumDigits, curIdx, curFigFileName));
        
        export_fig('-r172', strcat(curFigFileName, '.png'), curFig);
        export_fig('-r172', strcat(curFigFileName, '.eps'), curFig);
    end
end

function plotData(binEdges, synDists, ax, synClasses, classNames, varargin)
    opt = struct;
    opt.report = false;
    opt.classColors = num2cell( ...
        get(groot, 'defaultAxesColorOrder'), 2);
    opt = Util.modifyStruct(opt, varargin{:});
    
    %% Reporting
    if opt.report
        report = table;
        report.name = [{'Rest'}; classNames(:)];
        report.medianDistanceToSoma = accumarray( ...
            1 + synClasses, synDists, [], @median);
        report = report([2:end, 1], :);
        disp(report);
    end
    
    %% Plotting
    data = accumarray( ...
       [discretize(synDists, binEdges), 1 + synClasses], 1, ...
       [numel(binEdges) - 1, numel(classNames) + 1]);
    
    ax.TickDir = 'out';
    hold(ax, 'on');
    
    for curIdx = 1:numel(classNames)
        plotHist(ax, ...
            binEdges, data(:, 1 + curIdx), ...
            'EdgeColor', opt.classColors{curIdx});
    end
    
    leg = legend(ax, classNames);
    set(leg, 'Location', 'EastOutside', 'Box', 'off');
end

function plotHist(ax, edges, counts, varargin)
    histogram(ax, ...
        'BinEdges', edges, ...
        'BinCounts', counts, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1, ...
        varargin{:});
end
