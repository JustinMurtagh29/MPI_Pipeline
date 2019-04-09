% Author: 
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
% connectEM is run on the graph
% Hierarchical agglomeration is run on the graph
% Both are run here conservatively and the agglomerations from both 
% methods are here combined to get a larger, error-free agglomeration
% state

info = Util.runInfo();
borderSizeThr = 50;
segmentSizeThr = 100;
probThreshold = 0.99;
sizeThreshold = 1e6;
minScore = 0;
datasetName = 'H2_3_v2_U1_SubI';
% Load parameter
rootDir = ['/tmpscratch/sahilloo/data/' datasetName '/pipelineRun_mr2e_wsmrnet/'];
load(fullfile(rootDir,'allParameter.mat'))
Util.log('load data:')
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
graph = load([p.saveFolder 'graph.mat']);
maxSegId = Seg.Global.getMaxSegId(p);

% for skeletons
parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';

% output folder for saving new agglomeration state
folderName =  datestr(clock,30);
outputFolder = [p.saveFolder folderName '_agglomeration_NCHC/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

Util.log('cut and restrict graph based on NC thresholds')
graphCut = connectEM.cutGraphSimple(p, graph, segmentMeta, borderMeta, borderSizeThr, segmentSizeThr);
agglosNCEdges = graphCut.edges(graphCut.prob > probThreshold,:);
[~, agglosNC] = Graph.buildConnectedComponents(maxSegId, agglosNCEdges);
[agglosNC,agglosNCSizes] = doForAgglos(agglosNC, ['NC_' folderName], fullfile(outputFolder,'NC'),...
            graph.edges, segmentMeta, parameters, 100, true);

Util.log('do for HC graph')
graphHC = load(fullfile(rootDir, '29Nov2018_agglomeration/graph.mat'));
graphHC = graphHC.graph;
agglosHCEdges = graphHC.edges(graphHC.scores > minScore, :);
[~, agglosHC] = Graph.buildConnectedComponents(maxSegId, agglosHCEdges);
[agglosHC, agglosHCSizes] = doForAgglos(agglosHC, ['HC_' folderName], fullfile(outputFolder,'HC'),...
            graph.edges, segmentMeta, parameters, 100, true);

% NOTE (hack):
% very primitive,
% add: directionality,distance-threshold, typeEM information later
mergeEdges  = cat(1, agglosNCEdges, agglosHCEdges);
mergeEdges = unique(mergeEdges, 'rows','stable'); % get rid of duplicate edges
[~, agglos] = Graph.buildConnectedComponents(maxSegId, mergeEdges);
[agglos, agglosSize] = doForAgglos(agglos, ['NCHC_' folderName], outputFolder, graph.edges, segmentMeta, parameters, 100, false);

% Keep only agglomerates that have at least sizeThreshold million voxels
idx = agglosSize > sizeThreshold;
agglos = agglos(idx);
agglosSize = agglosSize(idx); 
clear idx

Util.log('save new agglomeration state')
Util.save(fullfile(outputFolder,'agglos.mat'),agglos, agglosSize, borderSizeThr,...
                    segmentSizeThr, probThreshold, sizeThreshold, minScore, info)

function [agglos, agglosSize] = doForAgglos(agglos, mappingName, outputFolder, edges, segmentMeta, parameters, count, saveFlag)
    if ~exist(outputFolder,'dir')
        mkdir(outputFolder)
    end
    agglosSize = cellfun(@(x) sum(segmentMeta.voxelCount(x)), agglos);
    [agglosSize, idx] = sort(agglosSize, 'descend');
    agglos = agglos(idx);
    clear idx
    Util.log('create mapping for new agglomeration state')
    WK.makeWKMapping(agglos, mappingName, outputFolder);
    
    Util.log('write skeletons for new agglomeration state')
    agglosOut = agglos(1:count);
    Superagglos.skeletonFromAgglo(edges, segmentMeta, ...
        agglosOut, 'agglos', outputFolder, parameters);
    if saveFlag
        Util.save(fullfile(outputFolder,[mappingName '_agglos.mat']),...
            agglos, agglosSize);
    end
end