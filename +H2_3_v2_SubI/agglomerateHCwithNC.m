% Author: 
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
% connectEM is run on the graph
% Hierarchical agglomeration is run on the graph
% Both are run here conservatively and the agglomerations from both 
% methods are here combined to get a larger, error-free agglomeration
% state

info = Util.runInfo();

% Load parameter
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/';
load(fullfile(rootDir,'allParameter.mat'))
Util.log('load data:')
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
graph = load([p.saveFolder 'graph.mat']);

% output folder for saving new agglomeration state
folderName =  datestr(clock,30);
outputFolder = [p.saveFolder folderName '_agglomeration_NCHC/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

Util.log('cut and restrict graph based on NC thresholds')
borderSizeThr = 50;
segmentSizeThr = 100;
probThreshold = 0.99;
sizeThreshold = 1e6;
graphCut = connectEM.cutGraphSimple(p, graph, segmentMeta, borderMeta, borderSizeThr, segmentSizeThr);
agglosNCEdges = graphCut.edges(graphCut.prob > probThreshold,:);

Util.log('load HC graph')
minScore = 0;
graphHC = load(fullfile(rootDir, '29Nov2018_agglomeration/graph.mat'));
graphHC = graphHC.graph;
agglosHCEdges = graphHC.edge(graphHC.score > minScore, :);

Util.log('catenate edges from two graphs')
% NOTE (hack): very primitive,
% add (i) directionality, (ii) distance-threshold (iii)typeEM information later
maxSegId = Seg.Global.getMaxSegId(p);
mergeEdges  = cat(1, agglosNCEdges, agglosHCEdges);
mergeEdges = unique(mergeEdges, 'rows'); % get rid of duplicate edges
[~, agglos] = Graph.buildConnectedComponents(maxSegId, mergeEdges);
agglosSize = cellfun(@(x) sum(segmentMeta.voxelCount(x)), agglos);
[agglosSize, idx] = sort(agglosSize, 'descend');
agglos = agglos(idx);
%{% Keep only agglomerates that have at least sizeThreshold million voxels
idx = agglosSize > sizeThreshold;
agglos = agglos(idx);
agglosSize = agglosSize(idx); %}
clear idx;
keyboard
Util.log('save new agglomeration state')
Util.save(fullfile(outputFolder,'agglos.mat'),agglos, agglosSize, borderSizeThr,...
                    segmentSizeThr, probThreshold, sizeThreshold, minScore, info)

Util.log('create mapping for new agglomeration state')
WK.makeWKMapping(agglos, ['NCHC:' folderName], outDir);

Util.log('write skeletons for new agglomeration state')
agglosOut = agglos(1:100);
parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';
tic;
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, ...
    agglosOut, 'agglos', outputFolder, parameters);
toc;

