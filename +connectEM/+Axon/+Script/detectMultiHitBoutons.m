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

[conn, syn] = connectEM.Connectome.load(param, connFile);

axons = load(conn.info.param.axonFile, 'axons');
axons = axons.axons;

[connDir, connName] = fileparts(connFile);
interSynFile = sprintf('%s_intersynapse.mat', connName);
interSynFile = fullfile(connDir, interSynFile);
interSyn = load(interSynFile);

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
curEdges = linspace(0, 10, 11);

curFig = figure();
curFig.Color = 'white';

curAx = axes(curFig);
hold(curAx, 'on');

curData = sortrows(annSynPairs, 'interSynDist');
curPrec = [1; cumsum(curData.sameBouton) ./ (1:height(curData))'; 0];
curRec = [0; cumsum(curData.sameBouton) / sum(curData.sameBouton); 1];

curAuc = curPrec .* curRec;
[~, curMaxAucIdx] = max(curAuc);
curData.interSynDist(curMaxAucIdx) / 1E3 %#ok

plot(curAx, curRec, curPrec);
curAx.XLim = [0, 1];
curAx.YLim = [0, 1];

xlabel(curAx, 'Recall');
ylabel(curAx, 'Precision');

%{
histogram(curAx, ...
    annSynPairs.interSynDist(annSynPairs.sameBouton) / 1E3, ...
    'BinEdges', curEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(curAx, ...
    annSynPairs.interSynDist(~annSynPairs.sameBouton) / 1E3, ...
    'BinEdges', curEdges, 'DisplayStyle', 'stairs', 'LineWidth', 2);

curAx.XLim = curEdges([1, end]);

xlabel(curAx, 'Intersynapse distance (µm)');
ylabel(curAx, 'Occurences');

curLeg = legend(curAx, 'Same bouton', 'Different bouton');
curLeg.Box = 'off';
%}

curAx.TickDir = 'out';

title(curAx, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);
