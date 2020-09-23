% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSyn = 0;

% Use WKW segmentation for speed
inParam = struct;
inParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
inParam.backend = 'wkwrap';

outParam = struct;
outParam.root = '/tmpscratch/amotta/l4/2019-11-01-connectome-segmentation/wkw/1';
outParam.backend = 'wkwrap';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
param.seg = inParam;

maxSegId = Seg.Global.getMaxSegId(param);
conn = connectEM.Connectome.load(param, connFile);

%% Build mapping
axons = conn.axons(conn.axonMeta.synCount >= minSyn);
dends = conn.dendrites(conn.denMeta.synCount >= minSyn);

confSegIds = intersect( ...
    cell2mat(axons), cell2mat(dends));
filter = @(agglos) cellfun(@(segIds) ...
    reshape(setdiff(segIds, confSegIds), [], 1), ...
    agglos, 'UniformOutput', false);

axons = filter(axons);
dends = filter(dends);

agglos = cat(1, axons, dends);
mappingIds = cellfun(@(ids) max([0, min(ids)]), agglos);

axonMappedIds = mappingIds(1:numel(axons));
dendMappedIds = mappingIds(numel(axons) + (1:numel(dends)));

assert(all(axonMappedIds(not(cellfun(@isempty, axons)))));
assert(all(dendMappedIds(not(cellfun(@isempty, dends)))));

mask = not(cellfun(@isempty, agglos));
mappingIds = mappingIds(mask);
agglos = agglos(mask);

mapping = Agglo.buildLUT(maxSegId, agglos, mappingIds);
mapping(not(mapping)) = find(not(mapping));
mapping = uint32(mapping);

%% Write output
outFile = fileparts(fileparts(outParam.root));
outFile = fullfile(outFile, 'meta.mat');

out = struct;
out.info = info;
out.axonMappedIds = axonMappedIds(:);
out.dendMappedIds = dendMappedIds(:);

Util.saveStruct(outFile, out);
Util.protect(outFile);

%% Build segmentation
wkwInit('new', outParam.root, 32, 32, 'uint32', 1);
Seg.Global.applyMappingToSegmentation(param, mapping, outParam);

% Create resolution pyramid
box = [1, 1, 1; (ceil(param.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(outParam, box, [], true);

% Compress segmentation
compressSegmentation(outParam);
