%% check for overlap in axon and dendrite agglos
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% get the overlapping eqivalence classes
m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
axonsPostQuery = m.axonsPostQuery;
m = load('E:\workspace\data\backed\L4\FocusEM\dendritesPostQuery.mat');
dendritesPostQuery = m.dendritesPostQuery;
combinedAgglos = {axonsPostQuery; dendritesPostQuery};
combinedAgglos = vertcat(combinedAgglos{:});
idx_ax = length(axonsPostQuery);
[ovAgglos, segId] = L4.Agglo.aggloOverlaps(combinedAgglos);

% get numbers on the overlap
ovAggloIds = unique(cell2mat(ovAgglos));
numOverlapping = length(ovAgglos);
numOV_axons = sum(ovAggloIds <= length(axonsPostQuery));
numOV_dendrites = sum(ovAggloIds > length(axonsPostQuery));
fraction_axons = numOV_axons/length(axonsPostQuery);
fraction_dendr = numOV_dendrites/length(dendritesPostQuery);

% overlap with multiple segments
l = cellfun(@length, segId);
minOverlap = 2;
ovAggloIdsM = ovAgglos(l >= minOverlap);
segIdM = segId(l >= minOverlap);
ovAggloIdsM = unique(cell2mat(ovAggloIdsM));
numOverlappingM = length(ovAggloIdsM);
numOV_axons_M = sum(ovAggloIdsM <= length(axonsPostQuery));
numOV_dendrites_M = sum(ovAggloIdsM > length(axonsPostQuery));
fraction_axons_M = numOV_axons_M/length(axonsPostQuery);
fraction_dendr_M = numOV_dendrites_M/length(dendritesPostQuery);

%% same for agglos above 5um only

% get the overlapping eqivalence classes
m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
m.isAbove5um(m.toIgnoreIdx) = false;
axonsPostQuery = m.axonsPostQuery(m.isAbove5um);
m = load('E:\workspace\data\backed\L4\FocusEM\dendritesPostQuery.mat');
dendritesPostQuery = m.dendritesPostQuery(m.isAbove5um);
combinedAgglos = {axonsPostQuery; dendritesPostQuery};
combinedAgglos = vertcat(combinedAgglos{:});
idx_ax = length(axonsPostQuery);
[ovAgglos, segId] = L4.Agglo.aggloOverlaps(combinedAgglos);

% get numbers on the overlap
ovAggloIds = unique(cell2mat(ovAgglos));
numOverlapping = length(ovAgglos);
numOV_axons = sum(ovAggloIds <= length(axonsPostQuery));
numOV_dendrites = sum(ovAggloIds > length(axonsPostQuery));
fraction_axons = numOV_axons/length(axonsPostQuery);
fraction_dendr = numOV_dendrites/length(dendritesPostQuery);

% overlap with multiple segments
l = cellfun(@length, segId);
minOverlap = 2;
ovAggloIdsM = ovAgglos(l >= minOverlap);
segIdM = segId(l >= minOverlap);
ovAggloIdsM = unique(cell2mat(ovAggloIdsM));
numOverlappingM = length(ovAggloIdsM);
numOV_axons_M = sum(ovAggloIdsM <= length(axonsPostQuery));
numOV_dendrites_M = sum(ovAggloIdsM > length(axonsPostQuery));
fraction_axons_M = numOV_axons_M/length(axonsPostQuery);
fraction_dendr_M = numOV_dendrites_M/length(dendritesPostQuery);

%% overlaps to nml

m = load('E:\workspace\data\backed\20170217_ROI\segmentMeta.mat', 'point');
points = m.point';
ax_off = length(axonsPostQuery);

% get random agglos with multiple overlap
ovAggloM = ovAgglos(l >= 2);
rng(1972017);
randIds = randperm(numel(ovAggloM), 50);
selAgglos = ovAggloM(randIds);
selSegId = segIdM(randIds);

% get respective segments for axons and dendrites and write them to
% skeleton
randAxAgglo = cell(length(randIds), 1);
randDeAgglo = cell(length(randIds), 1);
for i = 1:length(randIds)
    randAxAgglo{i} = axonsPostQuery{selAgglos{i}(selAgglos{i} <= ax_off)};
    randDeAgglo{i} = dendritesPostQuery{...
        selAgglos{i}(selAgglos{i} > ax_off) - ax_off};
end

randAgglos = [randAxAgglo, randDeAgglo]';
randAgglos = randAgglos(:);

% to skeleton
skel = Skeleton.fromMST( ...
    cellfun(@(ids) {points(ids, :)}, randAgglos(1:50)), [11.24, 11.24, 28]);
for i = 1:skel.numTrees()/2
    skel.names{2*i - 1} = sprintf('Pair%03d_Axon', i);
    skel.names{2*i} = sprintf('Pair%03d_Dendrite', i);
end
skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');

% add comments to overlapping nodes
for i = 1:skel.numTrees()/2
    [skel.nodesAsStruct{2*i - 1}(ismember(randAgglos{2*i - 1}, selSegId{i})).comment] = deal('overlapping');
%     [skel.nodesAsStruct{2*i}(ismember(randAgglos{2*i}, selSegId{i})).comment] = deal('overlapping');
end

%% size distribution of overlapping agglos
m = load('E:\workspace\data\backed\20170217_ROI\segmentMeta.mat', ...
    'voxelCount');
voxelCount = m.voxelCount;
ax_off = length(axonsPostQuery);

sizAD = zeros(length(ovAgglos), 2);
for i = 1:length(ovAgglos)
    ids = ovAgglos{i};
    sizAD(i, 1) = sum(voxelCount(axonsPostQuery{ ...
        ovAgglos{i}(ovAgglos{i} <= ax_off)}));
    sizAD(i, 2) = sum(voxelCount(dendritesPostQuery{ ...
        ovAgglos{i}(ovAgglos{i} > ax_off) - ax_off}));
end

scatter(sizAD(:, 1), sizAD(:,2), 'x');
xlabel('Axon agglo size (# voxels)')
ylabel('Dendrite agglo size (# voxels)')
a = gca;
a.TickDir = 'out';
a.FontSize = 24;
a.FontName = 'Arial';
title('Overlapping axon/dendrite pair sizes')
mA = max(sizAD(:,1));
mD = max(sizAD(:,2));
m = min(mA, mD);
line([0 m], [0 m]);

% zoom on
a.XLim = [0 1e7];
a.YLim = [0 1e7];

%% look at large dendrites

m = load('E:\workspace\data\backed\20170217_ROI\segmentMeta.mat', 'point');
points = m.point';
l = cellfun(@length, dendritesPostQuery);
size_t_upper = 6000;
size_t_lower = 4000;
largeDen = dendritesPostQuery(l > size_t_lower & l < size_t_upper);
skel = Skeleton.fromMST( ...
    cellfun(@(ids) {points(ids, :)}, largeDen), [11.24, 11.24, 28]); 
skel = Skeleton.setParams4Dataset(skel, 'ex145_ROI2017');
skel.write(sprintf('Dendrite_agglos_sizeRange_%d_%d.nml', ...
    size_t_lower, size_t_upper));