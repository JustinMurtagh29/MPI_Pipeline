% Detection of multi-hit / -synaptic boutons.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

synPairCount = 100;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

synT = connectEM.Connectome.buildSynapseTable(conn, syn);
synT.type = syn.synapses.type(synT.id);

axons = load(conn.info.param.axonFile, 'axons');
axons = axons.axons;

[connDir, connName] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse_v2.mat', connName);
interSynFile = fullfile(connDir, interSynFile);
interSyn = load(interSynFile);

boutonMetaFile = sprintf('%s_axonalBoutons_v1.mat', connName);
boutonMetaFile = fullfile(connDir, boutonMetaFile);
boutonMeta = load(boutonMetaFile);

%% Training data
% This section select N random synapses, and exports it together with one
% of its direct neighbours. I'll then decide whether or not a pair of
% synapses was originating from the same bouton.
interSyn.synCounts = cellfun(@numel, interSyn.synIds);
randAxonIdx = find(interSyn.synCounts > 1);

rng(0);
randAxonIdx = datasample( ...
    randAxonIdx, synPairCount, ...
    'weights', interSyn.synCounts(randAxonIdx));
randAxonIds = interSyn.axonIds(randAxonIdx);

randSynIds = nan(synPairCount, 2);
for curIdx = 1:synPairCount
    curAxonIdx = randAxonIdx(curIdx);
    
    curSynIds = interSyn.synIds{curAxonIdx};
    curSynToSynDists = interSyn.synToSynDists{curAxonIdx};
    curSynToSynDists(isinf(curSynToSynDists)) = 0;
    
    % Linearize axon
   [~, curRefIdx] = max(curSynToSynDists(:));
   [~, curRefIdx] = ind2sub(size(curSynToSynDists), curRefIdx);
    curRefDists = curSynToSynDists(:, curRefIdx);
    
    % Pick random synapse
    curSynIdx = randi(numel(curSynIds));
    curSynDists = curSynToSynDists(:, curSynIdx);
    
    % Pick random neighbour
    curLeftIdx = setdiff(find( ...
        curRefDists <= curRefDists(curSynIdx)), curSynIdx);
    curRightIdx = setdiff(find( ...
        curRefDists >= curRefDists(curSynIdx)), curSynIdx);
    
   [~, curTempIdx] = min(curSynDists(curLeftIdx));
    curLeftIdx = curLeftIdx(curTempIdx);
    
   [~, curTempIdx] = min(curSynDists(curRightIdx));
    curRightIdx = curRightIdx(curTempIdx);
    
    curNeighIdx = [curLeftIdx, curRightIdx];
    curNeighIdx = curNeighIdx(randi(numel(curNeighIdx)));
    
    randSynIds(curIdx, :) = curSynIds([curSynIdx, curNeighIdx]);
end

%% Generate skeleton for annotations
skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:synPairCount
    curAxonId = randAxonIds(curIdx);
    curParentId = conn.axonMeta.parentId(curAxonId);
    curTreeIds = skel.numTrees() + (1:3);
    
    skel = Superagglos.toSkel(axons(curParentId), skel);
    skel.names{end} = sprintf('Axon %d', curAxonId);
    
    curSynIds = randSynIds(curIdx, :);
    curSynAgglos = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    curSynAgglos = cellfun( ...
        @(ids) segPoints(ids, :), ...
        curSynAgglos, 'UniformOutput', false);
    
    skel = Skeleton.fromMST( ...
        curSynAgglos, param.raw.voxelSize, skel);
    skel.names((end - 1):end) = arrayfun( ...
        @(id) sprintf('Synapse %d', id), ...
        curSynIds, 'UniformOutput', false);
    
   [skel, curGroupId] = skel.addGroup(sprintf( ...
        'Axon %d, Synapses %d and %d', ...
        curAxonId, curSynIds(1), curSynIds(2)));
    skel = skel.addTreesToGroup(curTreeIds, curGroupId);
end

%% Loading annotations
annNmlFile = connectEM.Axon.Data.getDir('multiHitBoutonCalibration.nml');
annNmlFile = skeleton(annNmlFile);

% Process annotations
annSynPairs = annNmlFile.groups(:, {'id', 'name'});

curRegEx = '^axon (\d+), synapses (\d+) and (\d+) \((\w+)\)$';
annSynPairData = regexpi(annSynPairs.name, curRegEx, 'tokens', 'once');

annSynPairs(cellfun(@isempty, annSynPairData), :) = [];
annSynPairData(cellfun(@isempty, annSynPairData), :) = [];
annSynPairData = vertcat(annSynPairData{:});

annSynPairs.axonId = cellfun(@str2double, annSynPairData(:, 1));
annSynPairs.synIds = cellfun(@str2double, annSynPairData(:, 2:3));
annSynPairs.sameBouton = strcmpi(annSynPairData(:, 4), 'yes');
annSynPairs(:, {'name'}) = [];

% Find corresponding intersynapse distance
[~, annSynPairs.interSynIdx] = ismember( ...
    annSynPairs.axonId, interSyn.axonIds);
assert(all(annSynPairs.interSynIdx));

annSynPairs.interSynDist = nan(height(annSynPairs), 1);
for curIdx = 1:height(annSynPairs)
    curInterSynIdx = annSynPairs.interSynIdx(curIdx);
    curInterSynDists = interSyn.synToSynDists{curInterSynIdx};
    
    curDist = annSynPairs.synIds(curIdx, :);
   [~, curDist] = ismember(curDist, interSyn.synIds{curInterSynIdx});
    curDist = curInterSynDists(curDist(1), curDist(2));
    
    annSynPairs.interSynDist(curIdx) = curDist;
end

%% Plot results
clear cur*;
curPlot = false;
curEdges = linspace(0, 10, 11);

curData = sortrows(annSynPairs, 'interSynDist');
curPrec = [1; cumsum(curData.sameBouton) ./ (1:height(curData))'; 0];
curRec = [0; cumsum(curData.sameBouton) / sum(curData.sameBouton); 1];

curAuc = curPrec .* curRec;
[~, curMaxAucIdx] = max(curAuc);

% Determine optimal cut-off
optimalCutoff = curData.interSynDist(curMaxAucIdx);

if curPlot
    curFig = figure();
    curFig.Color = 'white';

    curAx = axes(curFig);
    hold(curAx, 'on');

    plot(curAx, curRec, curPrec);
    curAx.XLim = [0, 1];
    curAx.YLim = [0, 1];

    xlabel(curAx, 'Recall');
    ylabel(curAx, 'Precision');

    curAx.TickDir = 'out';

    title(curAx, ...
        {info.filename; info.git_repos{1}.hash}, ...
        'FontWeight', 'normal', 'FontSize', 10);
end

%% Define boutons
synIds = connectEM.Axon.getSynapses(conn, syn);

synIsSpine = cellfun( ...
    @(ids) syn.synapses.type(ids) == 'PrimarySpine', ...
    synIds, 'UniformOutput', false);

boutonIds = ...
    connectEM.Axon.clusterSynapsesIntoBoutons( ...
        synIds, interSyn, 'cutoffDist', optimalCutoff);

%% Plot average number of synapses per bouton
clear cur*;
curBinEdges = linspace(1, 3, 21);
curAxonClasses = axonClasses([1, 3]);

curMeanSynCounts = cellfun(@(boutonIds, isSpine) ...
    mean(accumarray(boutonIds, isSpine)), ...
    boutonIds, synIsSpine);

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [675, 375];

curAx = axes(curFig);
hold(curAx, 'on');

arrayfun( ...
    @(axClass) fixedHistogram(curAx, ...
        curMeanSynCounts(axClass.axonIds), ...
        'BinEdges', curBinEdges), ...
	curAxonClasses);

plot(curAx, ...
    repelem(1.2, 2), ylim(curAx), ...
    'Color', 'black', 'LineWidth', 2, 'LineStyle', '--');
plot(curAx, ...
    repelem(2.1, 2), ylim(curAx), ...
    'Color', 'black', 'LineWidth', 2, 'LineStyle', '--');

curVals = [1.2, 2.1];
curAx.XTick = union(curAx.XTick, curVals);

[~, curIds] = ismember(curVals, curAx.XTick);
curAx.XTickLabel(curIds) = {'CC', 'TC'};

set( ...
    curAx, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'XLim', curBinEdges([1, end]));

ylabel(curAx, 'Axons');
xlabel(curAx, 'Average number of spine synapses per bouton');

curLeg = legend(curAx, ...
    {curAxonClasses.title}, ...
    'Location', 'NorthEast');
curLeg.Box = 'off';

title( ...
    curAx, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot fraction of multihit boutons
clear cur*;
curBinEdges = linspace(0, 1, 21);
curAxonClasses = axonClasses([1, 3]);

curMultiHitFracs = cellfun(@(boutonIds, isSpine) ...
    mean(accumarray(boutonIds, isSpine) > 1), ...
    boutonIds, synIsSpine);

curFig = figure();
curFig.Color = 'white';
curFig.Position(3:4) = [675, 375];

curAx = axes(curFig);
hold(curAx, 'on');

arrayfun( ...
    @(axClass) fixedHistogram(curAx, ...
        curMultiHitFracs(axClass.axonIds), ...
        'BinEdges', curBinEdges), ...
	curAxonClasses);set( ...
    curAx, ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'XLim', curBinEdges([1, end]));

plot(curAx, ...
    repelem(0.19, 2), ylim(curAx), ...
    'Color', 'black', 'LineWidth', 2, 'LineStyle', '--');
plot(curAx, ...
    repelem(0.67, 2), ylim(curAx), ...
    'Color', 'black', 'LineWidth', 2, 'LineStyle', '--');

curVals = [0.19, 0.67];
curAx.XTick = union(curAx.XTick, curVals);

[~, curIds] = ismember(curVals, curAx.XTick);
curAx.XTickLabel(curIds) = {'CC', 'TC'};

ylabel(curAx, 'Axons');
xlabel(curAx, 'Fraction of boutons with multiple spine synapses');

curLeg = legend(curAx, ...
    {curAxonClasses.title}, ...
    'Location', 'NorthEast');
curLeg.Box = 'off';

title( ...
    curAx, {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Plot different identification criteria against each other
clear cur*;
curAxonIds = axonClasses(1).axonIds;

curConn = conn;
curConn.axonMeta.fullPriSpineSynDens = ...
    curConn.axonMeta.fullPriSpineSynCount ...
 ./ curConn.axonMeta.pathLen;

curMeanSynCounts = cellfun(@(boutonIds, isSpine) ...
    mean(accumarray(boutonIds, isSpine)), ...
    boutonIds, synIsSpine);
curMultiHitFracs = cellfun(@(boutonIds, isSpine) ...
    mean(accumarray(boutonIds, isSpine) > 1), ...
    boutonIds, synIsSpine);

curMedianBoutonVols = cellfun(@median, boutonMeta.boutonVols);
curMedianBoutonVols = curMedianBoutonVols / (1E3) ^ 3;

curGroupId = zeros(size(curAxonIds));
curGroupId(curMedianBoutonVols(curAxonIds) > 0.45) = 3;
curGroupId(curMeanSynCounts(curAxonIds) > 1.6) = 2;
curGroupId(curConn.axonMeta.fullPriSpineSynDens(curAxonIds) > 0.19) = 1;

curFig = figure();
curFig.Color = 'white';

curAx = subplot(1, 3, 1);
axis(curAx, 'square');
hold(curAx, 'on');

gscatter( ...
    curConn.axonMeta.fullPriSpineSynDens(curAxonIds), ...
    curMeanSynCounts(curAxonIds), curGroupId, [], [], [], 'off');
ylabel(curAx, 'Average spine synapses per bouton');

curAx = subplot(1, 3, 2);
axis(curAx, 'square');
hold(curAx, 'on');

gscatter( ...
    curConn.axonMeta.fullPriSpineSynDens(curAxonIds), ...
    curMultiHitFracs(curAxonIds), curGroupId, [], [], [], 'off');
ylabel(curAx, 'Fraction of multi-synaptic boutons');

curAx = subplot(1, 3, 3);
axis(curAx, 'square');
hold(curAx, 'on');

gscatter( ...
    curConn.axonMeta.fullPriSpineSynDens(curAxonIds), ...
    curMedianBoutonVols(curAxonIds), curGroupId, [], [], [], 'off');
ylabel(curAx, 'Median bouton volume (µm³)');

curAxes = flip(curFig.Children);

arrayfun(@(a) set(a, ...
    'TickDir', 'out', ...
    'XLim', [0, a.XLim(2)], ...
    'YLim', [0, a.YLim(2)]), curAxes);
xlabel(curAxes(1), 'Spine synapse density (µm^{-1})');

annotation( ...
    'textbox', [0, 0.9, 1, 0.1], ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'String', {info.filename; info.git_repos{1}.hash});

%% Look at discrepancies
curOutputDir = '/home/amotta/Desktop/tc-candidates';

rng(0);
curAxonIds = curAxonIds(curGroupId > 0);
curAxonIds = curAxonIds(randperm(numel(curAxonIds)));
curAxonIds = curAxonIds(1:10);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

mkdir(curOutputDir);

curNumDigits = ceil(log10(1 + numel(curAxonIds)));

for curIdx = 1:numel(curAxonIds)
    curId = conn.axonMeta.id(curAxonIds(curIdx));
    
    curSynIds = synIds{curId};
    curBoutonIds = boutonIds{curId};
    
    curAxon = conn.axons(curId);
    curSynapses = cellfun( ...
        @vertcat, ...
        syn.synapses.presynId(curSynIds), ...
        syn.synapses.postsynId(curSynIds), ...
        'UniformOutput', false);
    
    curNodes = [curAxon; curSynapses];
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curNodes, 'UniformOutput', false);
    
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, skel);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names(2:end) = arrayfun( ...
        @(id, boutonId) sprintf( ...
            'Synapse %d (Bouton %d)', id, boutonId), ...
        curSynIds, curBoutonIds, 'UniformOutput', false);
    curSkel.colors(2:end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', curNumDigits, curIdx, curId);
    curSkel.write(fullfile(curOutputDir, curSkelName));
end

%% Utilities
function varargout = fixedHistogram(ax, data, varargin)
    varargout = cell(1, nargout);
   [varargout{:}] = histogram( ...
        ax, data, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2, ...
        'FaceAlpha', 1, ...
        varargin{:});
end
