%% bouton agglomeration script
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo(false);
if ~ispc % assume on gaba
    pipelineFolder = ['/gaba/u/mberning/results/pipeline/' ...
        '20170217_ROI/'];
else % laptop
    pipelineFolder = 'E:\workspace\data\backed\20170217_ROI';
end

version = '_v3';
plotting = false;

%% load data

Util.log('Loading SVG data.');
stats = struct(); % stores result

m = load(fullfile(pipelineFolder, 'globalEdges.mat'), 'edges');
edges = m.edges;
m = load(fullfile(pipelineFolder, 'globalSynScores.mat'));
synScores = nan(size(edges));
synScores(m.edgeIdx, :) = m.synScores;
synT = -1.67;
stats.synT = synT;
mergeT = 0.98;
stats.mergeT = mergeT;
synEdgeIdx = max(synScores, [], 2) > -1.67;
m = load(fullfile(pipelineFolder, 'graph.mat'), 'prob', 'borderIdx');
prob = m.prob(~isnan(m.borderIdx));
m = load(fullfile(pipelineFolder, 'heuristicResult.mat'), 'myelinScore');
myelinScore = m.myelinScore;
m = load(fullfile(pipelineFolder, 'correspondencesNew2.mat'), 'corrEdges');
corrEdges = m.corrEdges;

% synapse predictions
synIdx = L4.Synapses.getSynapsePredictions(synScores, synT, true, prob, [], ...
    edges, myelinScore);

% restrict to edges at a 3 um margin to dataset border
m = load(fullfile(pipelineFolder, 'globalBorder.mat'), 'borderCoM');
borderCoM = m.borderCoM;
m = load(fullfile(pipelineFolder, 'allParameter.mat'), 'p');
p = m.p;
bbox = p.bbox;
margin = round(3000./[11.24; 11.24; 28]);
bboxM = bsxfun(@plus, bbox, [margin, -margin]);
isAtBorder = ~Util.isInBBox(borderCoM, bboxM);
synIdx(isAtBorder) = false;

%% get the boutons

Util.log('Get bouton agglomerates.');
sc = size(corrEdges, 1);

% separate bouton for each synapse
[bAgglo, postsynIdx, edgeIdx] = L4.Synapses.boutonAgglo( ...
    cat(1, edges, corrEdges), cat(1, synScores, nan(sc, 2)), ...
    cat(1, synIdx, false(sc, 1)), cat(1, prob, ones(sc, 1)), mergeT, 2);
boutons = table(bAgglo, postsynIdx, edgeIdx, 'VariableNames', ...
    {'agglo', 'postsynId', 'edgeIdx'});

% combine overlapping boutons
[bAgglo_nOv, postsynIdx_nOv, edgeIdx_nOv] = L4.Synapses.boutonAgglo(...
    cat(1, edges, corrEdges), cat(1, synScores, nan(sc, 2)), ...
    cat(1, synIdx, false(sc, 1)), cat(1, prob, ones(sc, 1)), mergeT, 2, true);
boutons_nOv = table(bAgglo_nOv, postsynIdx_nOv, edgeIdx_nOv, 'VariableNames', ...
    {'agglo', 'postsynId', 'edgeIdx'});

%% inspect 10 random bouton agglos

rng('shuffle');
numSamples = 10;
idx = randi(size(boutons_nOv, 1), numSamples, 1);
idx = [126940;33194;108442;16684;80401;66727;33644;14165;164172;86291];
m = load(fullfile(pipelineFolder, 'segmentMeta.mat'), 'point');
point = m.point';
skel = L4.Agglo.agglo2Nml(bAgglo_nOv(idx), point);
skel.names = arrayfun(@(x)sprintf('BoutonAgglo_%02d', x), 1:length(idx), ...
    'uni', 0);
for i = 1:numSamples
    connect_to = size(skel.nodes{i}, 1);
    for j = 1:length(postsynIdx_nOv{idx(i)})
        skel = skel.addNode(i, borderCoM(edgeIdx_nOv{idx(i)}(j),:), ...
            connect_to, [], 'synapse');
        skel = skel.addNode(i, point(postsynIdx_nOv{idx(i)}(j),:) , ...
            [], [], 'postsyn');
    end
end
% skel.write(['RandomBoutonAgglos' version '.nml'])
clear m

%% check for overlap of boutons with >5um axon agglos

Util.log('Calculating overlap statistics for axon agglos.');

% get axons >5um
if ~ispc
    % note that this does not contain the length -> AM
    m = load('/tmpscratch/mberning/axonQueryResults/postQueryAgglomerates.mat');
else
    m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
end
% m.isAbove5um(m.toIgnoreIdx) = false;
axonsPostQuery = m.axonsPostQuery(m.isAbove5um);
axonIds = cell2mat(axonsPostQuery);

m = max(cellfun(@max, bAgglo));
isInAxons = false(max([axonIds; m]), 1);
isInAxons(axonIds) = true;

boutonIsInAxons = cellfun(@(x)any(isInAxons(x)), bAgglo);
boutonOVIsInAxons = cellfun(@(x)any(isInAxons(x)), bAgglo_nOv);

stats.fraction_5um = sum(boutonIsInAxons)/length(boutonIsInAxons);
stats.total_5um = sum(boutonIsInAxons);
stats.fraction_5um_nOv = sum(boutonOVIsInAxons)/length(boutonOVIsInAxons);
stats.total_5um_nOv = sum(boutonOVIsInAxons);

%% check for overlap of boutons with all axon agglos

% get axons
if ~ispc
    m = load('/tmpscratch/mberning/axonQueryResults/postQueryAgglomerates.mat');
else
    m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
end

axonsPostQuery = m.axonsPostQuery;
axonIds = cell2mat(axonsPostQuery);

m = max(cellfun(@max, bAgglo));
isInAxons = false(max([axonIds; m]), 1);
isInAxons(axonIds) = true;

boutonIsInAxons = cellfun(@(x)any(isInAxons(x)), bAgglo);
boutonOVIsInAxons = cellfun(@(x)any(isInAxons(x)), bAgglo_nOv);

stats.fraction = sum(boutonIsInAxons)/length(boutonIsInAxons);
stats.total = sum(boutonIsInAxons);
stats.fraction_nOv = sum(boutonOVIsInAxons)/length(boutonOVIsInAxons);
stats.total_nOv = sum(boutonOVIsInAxons);

%% check for overlap with >5um axons including flight path collected agglos

m = load('/tmpscratch/mberning/axonQueryResults/axonQueryAnalysisBS.mat');
axonsNew = m.axonsNew;

b2agglo = L4.Synapses.getCollectedBoutons(axonsNew, bAgglo);
b2agglo_ov = L4.Synapses.getCollectedBoutons(axonsNew, bAgglo_nOv);

stats.fraction_flight = sum(b2agglo > 0)/length(b2agglo);
stats.total_flight = sum(b2agglo > 0);
stats.fraction_flight_nOv = sum(b2agglo_ov > 0)/length(b2agglo_ov);
stats.total_flight_nOv = sum(b2agglo_ov > 0);

%% display and save results
disp(stats)
if ispc
    outputFile = ['E:\workspace\data\backed\L4\Synapses\boutonAgglo' ...
        version '.mat'];
else
    outputFile = ['/u/bstaffle/data/L4/Synapses/boutonAgglos' ...
        version '.mat'];
end
if ~exist(outputFile, 'file')
    Util.log('Saving output to %s.', outputFile);
    save(outputFile, 'boutons', 'boutons_nOv', 'stats', 'info');
else
    warning('Output file %s already exists. Save the data manually.', ...
        outputFile);
end

%% unattached boutons to nml (from flight path collected agglos)

m = load(fullfile(pipelineFolder, 'segmentMeta.mat'), 'point');

numSamples = 20;
unattachedIdx = find(b2agglo_ov == 0);
rng('shuffle')
ridx = randi(length(unattachedIdx), numSamples, 1);
ridx = [50083; 3298; 68458; 77384; 106598; 3044; 43764; 93456; 39297; ...
    98361; 42780; 44005; 95815; 54635; 12446; 9251; 93192; 46799; 78650; 7913];
skel = L4.Agglo.agglo2Nml(bAgglo_nOv(unattachedIdx(ridx)), m.point');
skel.names = arrayfun(@(x)sprintf('MissedBoutonAgglo_%02d', x), ...
    1:length(ridx), 'uni', 0);
% skel.write(['MissedBoutonAgglos' version '.nml'])
clear m

%% bouton mapping

m = load(fullfile(pipelineFolder, 'globalSegSize.mat'));
maxSegId = length(m.segSize);

components = bAgglo_nOv(:);
toZero = setdiff(0:maxSegId, cell2mat(components));
components = cat(1, toZero, components);
% WK.makeWKMapping( components, ['boutons' version] );

%% bouton distance to dataset boundary

m = load(fullfile(pipelineFolder, 'segmentMeta.mat'), 'point');
point = m.point';
bComs = cell2mat(cellfun(@(x)mean(point(x,:), 1), bAgglo_nOv, 'uni', 0));
unattachedIdx = find(b2agglo_ov == 0);
bComs_unattached = cell2mat(cellfun(@(x)mean(point(x,:), 1), ...
    bAgglo_nOv(unattachedIdx), 'uni', 0));
d = Util.distToBoundary(bComs, bbox, [11.24, 11.24, 28]./1000);
d_unattached = Util.distToBoundary(bComs_unattached, bbox, ...
    [11.24, 11.24, 28]./1000);

if plotting
    % distance histogram (all boutons)
    figure;
    histogram(d);
    title('Bouton boundary distance')
    xlabel('Distance to boundary (\mu m)')
    ylabel('Bouton count (#)')
    Visualization.Figure.plotDefaultSettings()

    % distance histogram (missed boutons)
    figure;
    histogram(d_unattached);
    title('Bouton boundary distance')
    xlabel('Distance to boundary (\mu m)')
    ylabel('Bouton count (#)')
    Visualization.Figure.plotDefaultSettings()

    % cross-check: bouton center location in bbox
    bboxNM = bsxfun(@times, bbox, [11.24; 11.24; 28]./1000);
    bboxNM(:,2) = bboxNM(:,2) - bboxNM(:,1);
    bboxMNM = bsxfun(@times, bboxM, [11.24; 11.24; 28]./1000);
    bboxMNM(:,2) = bboxMNM(:,2) - bboxMNM(:,1);
    figure;
    linearize = @(x)x(:)';
    [~,ax, ~, h] = plotmatrix(bsxfun(@times, bComs_unattached, [11.24, 11.24, 28]./1000));
    rectangle(ax(1, 2), 'Position', linearize(bboxNM([2 1],:)));
    rectangle(ax(1, 3), 'Position', linearize(bboxNM([3 1],:)));
    rectangle(ax(2, 3), 'Position', linearize(bboxNM([3 2],:)));
    rectangle(ax(2, 1), 'Position', linearize(bboxNM([1 2],:)));
    rectangle(ax(3, 1), 'Position', linearize(bboxNM([1 3],:)));
    rectangle(ax(3, 2), 'Position', linearize(bboxNM([2 3],:)));
    rectangle(ax(1, 2), 'Position', linearize(bboxMNM([2 1],:)), 'EdgeColor', 'r');
    rectangle(ax(1, 3), 'Position', linearize(bboxMNM([3 1],:)), 'EdgeColor', 'r');
    rectangle(ax(2, 3), 'Position', linearize(bboxMNM([3 2],:)), 'EdgeColor', 'r');
    rectangle(ax(2, 1), 'Position', linearize(bboxMNM([1 2],:)), 'EdgeColor', 'r');
    rectangle(ax(3, 1), 'Position', linearize(bboxMNM([1 3],:)), 'EdgeColor', 'r');
    rectangle(ax(3, 2), 'Position', linearize(bboxMNM([2 3],:)), 'EdgeColor', 'r');
    title('Missed bouton locations')
    ax(1,1).YLabel.String = 'x (\mum)';
    ax(2,1).YLabel.String = 'y (\mum)';
    ax(3,1).YLabel.String = 'z (\mum)';
    ax(3,1).XLabel.String = 'x (\mum)';
    ax(3,2).XLabel.String = 'y (\mum)';
    ax(3,3).XLabel.String = 'z (\mum)';
    for i = 1:numel(ax)
        Visualization.Figure.plotDefaultSettings(ax(i));
    end
    Visualization.Figure.plotDefaultSettings(gca);
end

%% bouton synapse clustering examples

p = Gaba.getSegParameters('ex145_ROI2017');
m = load('/gaba/u/bstaffle/data/L4/Synapses/boutonAgglos_v3.mat');
boutons = m.boutons_nOv;
m = load(fullfile(p.saveFolder, 'globalEdges.mat'), 'edges');
edges = m.edges;
m = load(fullfile(p.saveFolder, 'globalNeuriteContinuityProb.mat'));
prob = m.prob;
m = load(fullfile(p.saveFolder, 'correspondencesNew.mat'));
edges = cat(1, edges, m.corrEdges);
prob(end+1:size(edges, 1)) = 1;
[synapses, boutonsNew] = L4.Synapses.clusterBoutonSynapsesLocalPostsynAgglo( ...
    boutons(1:1e4:end,:), edges, prob, 0.9, 2);

m = load(fullfile(p.saveFolder, 'globalBorder.mat'), 'borderCoM');
borderCom = m.borderCoM;
skel = L4.Synapses.boutonSynapses2Nml(synapses.syn_edgeIdx, borderCom);
skel.write('BoutonAggloPS.nml')

%% whole data bouton synapse clustering

[~, boutons] = L4.Synapses.clusterBoutonSynapsesLocalPostsynAgglo( ...
    boutons, edges, prob, 0.9, 2);

outputFile = ['/u/bstaffle/data/L4/Synapses/clusteredSynapses' ...
    version '.mat'];
if ~exist(outputFile, 'file')
    Util.log('Saving clustered synapses to %s.', outputFile);
    save(outputFile, 'boutons');
else
    warning('Output file %s already exists. Save the data manually.', ...
        outputFile);
end