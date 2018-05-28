% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSyn = 10;

% Use WKW segmentation for speed
inParam = struct;
inParam.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
inParam.backend = 'wkwrap';

outParam = struct;
outParam.root = '/tmpscratch/amotta/l4/2018-05-28-reconstruction-segmentation/wkw/1';
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
filter = @(agglos) cellfun( ...
    @(segIds) setdiff(segIds, confSegIds), ...
    agglos, 'UniformOutput', false);

axons = filter(axons);
axons(cellfun(@isempty, axons)) = [];

dends = filter(dends);
dends(cellfun(@isempty, dends)) = [];

confSegIds = num2cell(confSegIds(:), 2);
agglos = cat(1, axons, dends, confSegIds);

mappingIds = cellfun(@(ids) ids(1), agglos);
mapping = uint32(Agglo.buildLUT(maxSegId, agglos, mappingIds));

%% Build segmentation
wkwInit('new', outParam.root, 32, 32, 'uint32', 1);
Seg.Global.applyMappingToSegmentation(param, mapping, outParam);

% Create resolution pyramid
thisBBox = [1, 1, 1; (ceil(param.bbox(:, 2) ./ 1024) .* 1024)']';
createResolutionPyramid(outParam, thisBBox, [], true);

% Compress segmentation
compressSegmentation(outParam);

%% Write info
infoFile = fileparts(fileparts(outParam.root));
infoFile = fullfile(infoFile, 'info.mat');

Util.save(infoFile, info);
Util.protect(infoFile);
