% script to generate the dendrite gt whole cell gallery
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();

%% load data

p = Gaba.getSegParameters('ex145_ROI2017');
% p.agglo.dendriteAggloFile = fullfile(p.agglo.saveFolder, ...
%     'dendrites_03_v2_splitmerged.mat');
p.agglo.dendriteAggloFile = fullfile(p.agglo.saveFolder, ...
    'dendrites_04.mat');
m = load(p.svg.edgeFile);
edges = m.edges;
maxSegId = max(m.maxEdgeIdx);
m = load(p.agglo.dendriteAggloFile);
dendrites = m.dendrites;
m = load(p.agglo.axonAggloFile);
axons = m.axons;
m = load(p.svg.segmentMetaFile, 'point');
point = m.point';

% get skeletons from connectEM/evaluationData/wholeCellGT
thisFolder = fileparts(mfilename('fullpath'));
skelPath = fullfile(fileparts(fileparts(thisFolder)), 'evaluationData', ...
    'wholeCellGT');
skels = skeleton.loadSkelCollection(skelPath, [], true);

% find larger tree and assume it is the dendrite
treeSize = cellfun(@(skel)cell2mat(cellfun(@(trNodes)size(trNodes, 1), ...
    skel.nodes, 'uni', 0)), skels, 'uni', 0);
skels = cellfun(@(skel, tr)skel.deleteTrees(tr == min(tr)), ...
    skels, treeSize, 'uni', 0);

% get seed coords for exclusion of bbox around soma
seedCoords = zeros(length(skels), 3);
for i = 1:length(skels)
    [~, idx] = skels{i}.getNodesWithIDs(1);
    seedCoords(i, :) = skels{i}.nodes{1}(idx, 1:3);
end

% restrict to center
boundaryNM = 5e3; % um
bboxR = Util.addBorder(p.bbox, -ceil(boundaryNM./[11.24; 11.24; 28]));
skelsR = cellfun(@(x)x.restrictToBBox(bboxR, [], false), skels, 'uni', 0);

% delete parts around seed coord
l = 5e3; % um on each side of seed
bbox = ceil(l./p.raw.voxelSize(:));
bbox = [-bbox, bbox];
bboxExclSoma = cell(length(skels), 1);
for i = 1:length(skels)
    nodes = skelsR{i}.nodes{1}(:, 1:3);
    bboxExclSoma{i} = bsxfun(@plus, seedCoords(i, :)', bbox);
    toDel = Util.isInBBox(nodes, bboxExclSoma{i});
    skelsR{i} = skelsR{i}.deleteNodes(1, toDel);
end


%% get seg ids (needed in several other functions this I do it here)

segIds = Skeleton.getSegmentIdsOfSkelCellArray(p, skelsR);


%% get soma/nuclei and add them to agglos

Util.log('Adding soma agglos to dendrites.');
somaAgglos = load(fullfile(p.agglo.saveFolder,'center_somas.mat'));
somaAgglos = somaAgglos.result(:,2);
idx = ~cellfun(@isempty, somaAgglos);
denIds = Superagglos.getSegIds(dendrites);
denLUT = Agglo.buildLUT(maxSegId, denIds);
somaToDen = cellfun(@(x)denLUT(x), somaAgglos, 'uni', 0);
somaToDen = cellfun(@(x)x(x > 0), somaToDen, 'uni', 0);
somaOv = cellfun(@(x)tabulate(x), somaToDen, 'uni', 0);
somaOv(idx) = cellfun(@(x)x(x(:,2) > 0, :), somaOv(idx), 'uni', 0);
noDenSegs = cell(length(somaAgglos), 1);
noDenSegs(idx) = cellfun(@(y)arrayfun(@(x)size(x.nodes, 1), ...
    dendrites(y(:,1))), somaOv(idx), 'uni', 0);
somaOv(idx) = cellfun(@(x,y) cat(2, x, y), somaOv(idx), ...
    noDenSegs(idx), 'uni', 0);
fullySomaOverlapping = cell(length(somaAgglos), 1);
fullySomaOverlapping(idx) = cellfun(@(x)x(x(:,2) == x(:,4), 1), ...
    somaOv(idx), 'uni', 0);
idx = cellfun(@mode, somaToDen);
for i = 1:length(idx)
    if ~isnan(idx(i))
        % simply merge full soma into dendrite with most overlapping
        % segments
        dendrites(idx(i)).nodes = cat(1, dendrites(idx(i)).nodes, ...
            [point(somaAgglos{i}, :), somaAgglos{i}]);
        
%         % merge all fully soma overlapping dendrite agglos into most
%         % overlapping dendrite agglo
%         dendrites(idx(i)).nodes = cat(1, dendrites(idx(i)).nodes, ...
%             dendrites(fullySomaOverlapping{i}).nodes);
    end
end

% delete dendrite agglos that were merged into largest overlapping one
dendrites(cell2mat(fullySomaOverlapping)) = [];


%% make gallery

Util.log('Calculting agglo overlaps & metrics.');

% skeleton agglo overlap
skelToDendrites = L4.Agglo.aggloSkelOverlap(segIds, [], dendrites);
ovT = 3; % only agglos overlapping with >= ovT nodes are kept
skelToDendritesT = cellfun(@(x)x(x(:,2) >= ovT, :), skelToDendrites, ...
    'uni', 0);
[ovL_old, recL_old] = L4.Agglo.aggloSkelPathLengthOverlap(skelsR, segIds, ...
    cellfun(@(x)x(:,1), skelToDendritesT, 'uni', 0), dendrites);
[ovL, recL] = L4.Agglo.aggloSkelPathLengthOverlap2(skelsR, segIds, ...
    cellfun(@(x)x(:,1), skelToDendritesT, 'uni', 0), dendrites, ovT - 1);

% erl (note that if there is a node placed in a different segment along a
% path the erl stops there although the large agglo might just miss that
% single node)
eClasses = Superagglos.getSegIds(cat(1, axons, dendrites));
CRC = L4.Segmentation.ERL.calculateCRC(skelsR, segIds, edges, [], ...
    eClasses, true);
[erl, skel_erl] = L4.Segmentation.ERL.calculateERL(CRC);
erl = erl / 1000; % um
skel_erl = skel_erl ./ 1000; % um

% sort overlapping dendrite agglos by length
for i = 1:length(skelToDendritesT)
    [ovL{i}, sidx] = sort(ovL{i}, 'descend');
    skelToDendritesT{i} = skelToDendritesT{i}(sidx, :);
%     ovDenLength{i} = ovDenLength{i}(sidx);
end

% split length (excluding nodes outside the segmentation bbox)
L = cellfun(@(x)x.pathLength()./1000, skelsR);
skel_wavgL = arrayfun(@(x, l_tot) ...
    L4.Agglo.lengthWeightedAvg(x{1}, l_tot), ovL, L);
wavgL = L4.Agglo.lengthWeightedAvg(cell2mat(ovL), sum(L));

% splits based solely on the skelToAgglos mapping
numDenAgglos = cellfun(@(x)size(x, 1), skelToDendritesT);

% count nodes not in dendrite agglos
denIds = Superagglos.getSegIds(dendrites);
tmp = cell2mat(denIds);
unaccountedIds = cellfun(@(x)setdiff(x{1}, [tmp; 0; -1]), segIds, 'uni', 0);

% group unaccounted ids by axons and rest
axonsLUT = Agglo.buildLUT(maxSegId, Superagglos.getSegIds(axons));
axonAggloIdx = cellfun(@(x)axonsLUT(x), unaccountedIds, 'uni', 0);
skelToAxons = cellfun(@(x)setdiff(x, 0), axonAggloIdx, 'uni', 0);
numAxonAgglos = cellfun(@length, skelToAxons);
numRest = cellfun(@(x)sum(x == 0), axonAggloIdx);

% add both kinds of splits and calculate split distance
numSplits = numDenAgglos; % + numAxonAgglos + numRest - 1;
d_s = sum(L)./sum(numSplits);
skel_d_s = L./numSplits;

% results
result.d_s = d_s;
result.wavgL = wavgL;
result.L = sum(L);
result.skel = table(L, recL, numSplits, skel_d_s, skel_wavgL, ovL, ...
    'VariableNames', {'Length', 'RecoveredLength', 'NumSplits', ...
    'SplitLength', 'WeighAvgL', 'SingleAggloOverlapL'});


%% dendrite agglo overlap gallery

Util.log('Creating overlap gallery.');

outputFolder = '/gaba/u/bstaffle/data/L4/Dendrites/wholeCellAggloEval/gallery/';
noPlot = false;
ovSkels = L4.Agglo.skelOverlapGallery(outputFolder, skels, dendrites, ...
    skelToDendritesT, ovT, [], noPlot, L, ovL, zeros(length(skels), 1), ...
    skel_erl, skel_wavgL, recL);

save(fullfile(outputFolder, 'skelAggloOverlap.mat'), 'skels', 'skelsR', ...
    'segIds', 'ovL', 'skelToDendrites', 'L', 'd_s', ...
    'skel_d_s', 'ovSkels', 'info')


%% get nodes in agglo

% skelNo = 16;
% aggloNo = 1372274;
% aggloIds = Superagglos.getSegIds(dendrites(aggloNo));
% nodes = skelsR{skelNo}.nodes{1}(ismember(segIds{skelNo}{1},aggloIds{1}), 1:3);

%% ov skels to nml for corrections

% for i = 1:length(ovSkels)
%     ovSkels{i}.write(sprintf('Skel%02d_%s_denOverlaps.nml', i, ...
%         ovSkels{i}.filename));
% end

%% all overlaps for proofreading

% ovSkelsAll = L4.Agglo.aggloSkelOverlap2Nml(skelsR, segIds, dendrites, ...
%     axons, point);

%% merger detection
m = load(p.svg.borderMetaFile, 'borderCoM');
borderCom = m.borderCoM;
distT = 3e3;
mergers = L4.WholeCellGT.skelDistMergerDet( skelsR(1), segIds, dendrites, ...
    edges, borderCom, skelToDendritesT, distT, p.raw.voxelSize, bboxR, ...
    bboxExclSoma );

% write tree with merger borders to skeleton
idx = 1;
skel = skelsR{idx};
skel = skel.splitCC();
skel.colors(1:end) = {[1 0 0 1]};
nT = skel.numTrees();
skel.names(1:nT) = arrayfun(@(x)sprintf('GT_part%d', x), 1:nT, 'uni', 0);
skel = L4.Agglo.agglo2Nml(Superagglos.getSegIds( ...
    dendrites(skelToDendritesT{idx}(:,1))), point, skel);
skel.names(nT + 1:end) = arrayfun(@(x)sprintf('Agglo_%d', x), ...
    skelToDendritesT{idx}(:,1), 'uni', 0);
skel.colors(nT+1:end) = {[0 1 0 1]};
skel = skel.addTree([], mergers.borderCom{idx}{1});
skel.names{end} = 'PossibleLocations';
skel.colors{end} = [0 0 1 1];