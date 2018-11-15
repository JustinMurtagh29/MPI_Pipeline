% Manually inspect configurations where an excitatory axon makes contact
% with exactly two spine synapses of a dendrite.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
outputDir = '/home/amotta/Desktop/double-spine-contacts';

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

%% Find double spine contacts
clear cur*;
[axonShT, ~, graph.uniId] = unique(graph(:, {'shIdx', 'axonId'}), 'rows');
axonShT.dendId = shT.dendId(axonShT.shIdx);

[axonDendT, ~, curShIndices] = unique( ...
    axonShT(:, {'axonId', 'dendId'}), 'rows');
axonDendT.shInd = accumarray( ...
    curShIndices, axonShT.shIdx, [], @(ids) {ids});
axonDendT = axonDendT(cellfun(@numel, axonDendT.shInd) == 2, :);

%% Export random examples to webKnossos
clear cur*;
exportRange = 1:10;

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

        curAxon = conn.axons{curAxonId};
        curDend = conn.dendrites{curDendId};
        curShT = shT(curAxonDendT.shInd{curIdx}, :);

        curSkel = skeleton();

        curSkel = Skeleton.fromMST( ...
            points(curAxon, :), param.raw.voxelSize, curSkel);
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
