% use spine attached dend agglos to write segmentation on wK
% also add axons agglomerates to same segmentation

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;
p = param;

axonFile = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/aggloMat/20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results.mat';
dendriteFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results_auto-spines_v3.mat');

%% load data
Util.log('Loading data...')
segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
maxSegId = segmentMeta.maxSegId;

Util.log('Dendrite version: %s',dendriteFile)
dendrites = Util.load(dendriteFile, 'dendrites');
dendAgglos = arrayfun(@(x) SuperAgglo.toBundle(maxSegId,x), dendrites);
dendAgglos = arrayfun(@(x) x.agglomerate_to_segments, dendAgglos, 'uni',0);

Util.log('Axon version: %s',axonFile)
axons = Util.load(axonFile, 'axons');
axonAgglos = arrayfun(@(x) x.segIds, axons, 'uni',0);

%% segmentation to WKW
Util.log('Write new segmentation based on agglos')
segOut = struct;
segOut.root = fullfile(p.saveFolder, 'aggloState', ...
    ['20191227T134319_ga_20191224T235355optimParams' '_seg_axons_results_auto-spines_v3'], '1');

mkdir(segOut.root)
dataType = 'uint32';
wkwInit('new',segOut.root,32, 32, dataType, 1);
segOut.backend = 'wkwrap';

agglosSorted = vertcat(dendAgglos,axonAgglos);
mapping = connectEM.createLookup(segmentMeta, agglosSorted);
Seg.Global.applyMappingToSegmentation(p, mapping, segOut);
thisBBox = [1, 1, 1; (ceil(p.bbox(:, 2) ./ 1024) .* 1024)']';


createResolutionPyramid(segOut, thisBBox, [], true);

compressSegmentation(segOut)

