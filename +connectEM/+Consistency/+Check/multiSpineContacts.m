% Manually inspect configurations where an excitatory axon makes contact
% with at least two spine synapses of a dendrite.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

outputDir = '';
annotationDir = fullfile( ...
    fileparts(mfilename('fullpath')), ...
    'annotations', 'multi-spine-contacts');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
points = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);

% Loading spine head agglomerates
shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

% Loading augmented graph
graph = Graph.load(rootDir);
graph(~graph.borderIdx, :) = [];
graph(:, {'prob'}) = [];

borderAreas = fullfile(rootDir, 'globalBorder.mat');
borderAreas = load(borderAreas, 'borderArea2');
borderAreas = borderAreas.borderArea2;

graph.borderArea = borderAreas(graph.borderIdx);
clear borderAreas;

%% Assign spine heads to dendrites
shT = table;
shT.id = reshape(1:numel(shAgglos), [], 1);
shT.agglo = shAgglos;

% NOTE(amotta): The comment and code block below was reproduced from
% connectEM.Connectome.plotSynapseSizeConsistency

% NOTE(amotta): Find spine heads that overlap with axon agglomerates.
% This is a sign that something went wrong. We don't want to use these
% data points in our analysis.
%   By removing the spine head we make sure that the corresponding
% synapse is not detected as being onto a spine head, and thus won't
% make it into the axon-spine interface table.
curAxonLUT = Agglo.buildLUT(maxSegId, conn.axons);
shT = shT(cellfun(@(ids) ~any(curAxonLUT(ids)), shT.agglo), :);

% Find dendrite IDs for spine heads
curDendLUT = Agglo.buildLUT(maxSegId, conn.dendrites);

shT.dendId = cellfun( ...
    @(ids) setdiff(curDendLUT(ids), 0), ...
    shT.agglo, 'UniformOutput', false);

% Only consider spine heads that intersect with exactly one dendrite
shT = shT(cellfun(@isscalar, shT.dendId), :);
shT.dendId = cell2mat(shT.dendId);

%% Find spine heads touched by each axon
clear cur*;

curShLUT = Agglo.buildLUT(maxSegId, shT.agglo);
graph.shIdx = curShLUT(graph.edges);
graph = graph(xor(graph.shIdx(:, 1), graph.shIdx(:, 2)), :);
graph.shIdx = max(graph.shIdx, [], 2);

curAxonLUT = Agglo.buildLUT(maxSegId, conn.axons);
% NOTE(amotta): At this point we already know that exactly one of the two
% segments of the remaining edges is occupied by a spine head, and that the
% spine heads in question do not overlap with axons.
graph.axonId = max(curAxonLUT(graph.edges), [], 2);
graph = graph(graph.axonId > 0, :);

%% Find multi-spine contacts
clear cur*;
[axonShT, ~, graph.axonShId] = unique( ...
    graph(:, {'shIdx', 'axonId'}), 'rows');
axonShT.dendId = shT.dendId(axonShT.shIdx);
axonShT.area = accumarray(graph.axonShId, graph.borderArea);

[axonDendT, ~, curShIndices] = unique( ...
    axonShT(:, {'axonId', 'dendId'}), 'rows');
axonDendT.shInd = accumarray( ...
    curShIndices, axonShT.shIdx, [], @(ids) {ids});
axonDendT = axonDendT(cellfun(@numel, axonDendT.shInd) > 1, :);

%% Export random examples to webKnossos
% Note that this section is essentialy skipped if `outputDir` is not set
clear cur*;
exportRange = 1:100;

% HACK(amotta): Only run exports if output directory was set
if isempty(outputDir); exportRange = []; end

% NOTE(amotta): Only look at spine synapses from excitatory axons
curAxonDendT = axonDendT(ismember( ...
    axonDendT.axonId, axonClasses(1).axonIds), :);

% Randomize order (reproducibly)
rng(0);
curAxonDendT = curAxonDendT(randperm(height(curAxonDendT)), :);
curAxonDendT = curAxonDendT(exportRange, :);

curNumDigits = ceil(log10(1 + numel(exportRange)));

for curIdx = 1:numel(exportRange)
    try
        curAxonId = curAxonDendT.axonId(curIdx);
        curDendId = curAxonDendT.dendId(curIdx);

        curAxonSegIds = conn.axons{curAxonId};
        curDend = conn.dendrites{curDendId};
        curShT = shT(curAxonDendT.shInd{curIdx}, :);

        curSkel = skeleton();

        curSkel = Skeleton.fromMST( ...
            points(curAxonSegIds, :), param.raw.voxelSize, curSkel);
        curSkel.names{end} = sprintf('Axon %d', curAxonId);
        curSkel.colors{end} = [1, 0, 0, 1];

        curSkel = Skeleton.fromMST( ...
            points(curDend, :), param.raw.voxelSize, curSkel);
        curSkel.names{end} = sprintf('Dendrite %d', curDendId);
        curSkel.colors{end} = [0, 0, 1, 1];

        curSkel = Skeleton.fromMST( ...
            cellfun(@(ids) {points(ids, :)}, curShT.agglo), ...
            param.raw.voxelSize, curSkel);
        curSkel.names(3:end) = arrayfun(@(id) sprintf( ...
            'Spine head %d', id), curShT.id, 'UniformOutput', false);
        curSkel.colors(3:end) = {[0, 1, 0, 1]};

        curSkel = Skeleton.setParams4Pipeline(curSkel, param);
        curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
        
        curName = sprintf( ...
            '%d_axon-%d_dendrite-%d.nml', ...
            exportRange(curIdx), curAxonId, curDendId);
        curSkel.write(fullfile(outputDir, curName));
    catch err
        % NOTE(amotta): Try to prevent out-of-memory crashed.
        disp(err)
    end
end

%% Evaluate annotated contact sites
clear cur*;

curNmlFiles = dir(fullfile(annotationDir, '*.nml'));
curNmlFiles = {curNmlFiles(~[curNmlFiles.isdir]).name};

annT = table;
annT.nmlFile = reshape(curNmlFiles, [], 1);

for curIdx = 1:height(annT)
    curNmlName = annT.nmlFile{curIdx};
    curNmlPath = fullfile(annotationDir, curNmlName);
    
    % HACK(amotta): Extract axon and dendrite ID from file name
    curIds = '\d+_axon-(\d+)_dendrite-(\d+)\.nml';
    curIds = regexpi(curNmlName, curIds, 'tokens', 'once');
    
    curIds = cellfun(@str2double, curIds, 'UniformOutput', false);
   [curAxonId, curDendId] = deal(curIds{:});
   
    curExpectedShCount = sum( ...
        axonShT.axonId == curAxonId ...
      & axonShT.dendId == curDendId);
   
    curNml = slurpNml(curNmlPath);
    curTrees = NML.buildTreeTable(curNml);
    
    curAxonMask = startsWith(curTrees.name, 'Axon');
    curShMask = startsWith(curTrees.name, 'Spine head');
    
    assert(sum(curShMask) == curExpectedShCount);
    assert(all(curAxonMask | curShMask));
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes.coord = curNodes.coord + 1;
    
    % NOTE(amotta): Agglomerate one is the axon, followed by the spine head
    % agglomerates in the same order as in `curTrees`.
   [~, curNodes.aggloId] = ismember( ...
       curNodes.treeId, curTrees.id(curShMask));
    assert(all(ismember( ...
        curNodes.treeId(~curNodes.aggloId), ...
        curTrees.id(curAxonMask))));
    
    curNodes.aggloId = curNodes.aggloId + 1;
    curNodes.segId = Seg.Global.getSegIds(param, curNodes.coord);
    
    curAgglos = accumarray( ...
        curNodes.aggloId, curNodes.segId, [], @(segIds) {segIds});
    curLUT = Agglo.buildLUT(maxSegId, curAgglos);
    Agglo.check(curAgglos);
    
    % Calculate ASI areas
    curGraph = graph;
    curGraph.shId = curLUT(curGraph.edges);
    
    curGraph = curGraph( ...
        any(curGraph.shId == 1, 2) ...
      & any(curGraph.shId>= 2, 2), :);
    curGraph.shId = max(curGraph.shId, [], 2) - 1;
  
    curSh = regexpi( ...
        curTrees.name(curShMask), ...
        'spine head (\d+) \((\w+)\)', ...
        'tokens', 'once');
    
    curShT = vertcat(curSh{:});
    curShT = cell2table(curShT, 'VariableNames', {'shId', 'isSyn'});
    
    curShT.shId = cellfun(@str2double, curShT.shId);
    curShT.isSyn = strcmpi(curShT.isSyn, 'yes');
    
    curShT.area = accumarray( ...
        curGraph.shId, curGraph.borderArea, [height(curShT), 1]);
    annT.areas{curIdx} = curShT.area(curShT.isSyn);
end

%% Plot results
clear cur*;

curCount = cellfun(@numel, annT.areas);
curGroup = repelem(curCount, curCount);
curAreas = cell2mat(annT.areas);


% TODO(amotta): Find out what the source of these NaNs is and make sure it
% does not systematically affect the result.
curGroupMeanLogs = accumarray( ...
    curGroup, curAreas, [], @(v) mean(log10(v), 'omitnan'));

curMeanRatios = ...
    mean(curAreas(curGroup == 2), 'omitnan') ...
  / mean(curAreas(curGroup == 1), 'omitnan') %#ok

curMeanLogRatio = ...
    10^curGroupMeanLogs(2) ...
  / 10^curGroupMeanLogs(1) %#ok


% Histogram
curFig = figure;
curFig.Position(3:4) = [385, 335];
curFig.Color = 'white';
curAx = axes(curFig);
hold(curAx, 'on');

curBins = linspace(-2, 1, 7);
curHistOne = histogram(curAx, ...
    log10(curAreas(curGroup == 1)), ...
    'BinEdges', curBins, 'DisplayStyle', 'stairs', ...
    'LineWidth', 2, 'Normalization', 'probability');
curHistTwo = histogram(curAx, ...
    log10(curAreas(curGroup == 2)), ...
    'BinEdges', curBins, 'DisplayStyle', 'stairs', ...
    'Normalization', 'probability', 'LineWidth', 2);

curColors = get(groot, 'defaultAxesColorOrder');

plot(curAx, ...
    repelem(curGroupMeanLogs(1), 1, 2), curAx.YLim, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', curColors(1, :));

plot(curAx, ...
    repelem(curGroupMeanLogs(2), 1, 2), curAx.YLim, ...
    'LineWidth', 2, 'LineStyle', '--', 'Color', curColors(2, :));

xlabel(curAx, 'log10(axon-spine interface area (µm²))');
ylabel(curAx, 'Probability');

curAx.XLim = curBins([1, end]);
curAx.TickDir = 'out';
curLeg = legend(curAx, { ...
    '1 spine synapse per connection', ...
    '2 spine synapses per connection'}, ...
    'Location', 'SouthOutside');
curLeg.Box = 'off';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);


% Box plot
curFig = figure();
curFig.Position(3:4) = [290, 420];
curFig.Color = 'white';
curAx = axes(curFig);

boxplot(curAx, log10(curAreas), curGroup, 'width', 0.8)
xlabel(curAx, 'Spine synapses per connection');
ylabel(curAx, 'log10(axon-spine interface area (µm²))');

curAx.TickDir = 'out';
curAx.Box = 'off';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Look at synapse size similarity
clear cur*;

curCv = @(a) std(a) / mean(a);
curSaSdMask = cellfun(@numel, annT.areas) > 1;
curSaSdCvs = cellfun(curCv, annT.areas(curSaSdMask));

curRandCvs = cell2mat(annT.areas(~curSaSdMask));
curRandCvs = curRandCvs(randi(numel(curRandCvs), [1000, 2]));
curRandCvs = cellfun(curCv, num2cell(curRandCvs, 2));

% Plot
curFig = figure;
curFig.Color = 'white';
curFig.Position(3:4) = [450, 305];

curAx = axes(curFig);
curAx.TickDir = 'out';
hold(curAx, 'on');

curBins = linspace(0, 1.5, 6);
curHist = @(data) histogram( ...
    curAx, data, 'BinEdges', curBins, ...
    'DisplayStyle', 'stairs', 'LineWidth', 2, ...
    'Normalization', 'probability');

curHist(curSaSdCvs);
curHist(curRandCvs);

xlim(curAx, curBins([1, end]));
xlabel(curAx, 'Coefficient of variation');
ylabel(curAx, 'Probability');

curLeg = legend(curAx, { ...
    'Same-axon same-dendrite', ...
    'Random pairs of spine synapses'}, ...
    'Location', 'NorthEast');
curLeg.Box = 'off';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
