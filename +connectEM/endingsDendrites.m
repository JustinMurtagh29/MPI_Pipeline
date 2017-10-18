
% load('/gaba/scratch/mberning/20170404T172554_agglomeration/aggloNew.mat','dendritesFinal');
dendrites = load('/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_April.mat')
dendrites = dendrites.dendrites;
dendrites_v2 = load('/u/mberning/results/pipeline/20170217_ROI/aggloState/dendrites_03_v2.mat')
dendrites_v2 = dendrites_v2.dendrites(dendrites_v2.ind;
dendrites_v2 = dendrites_v2.

options.latentScore = 0.7;
options.segDirScore = 0.9;
options.neuriCScore = 0.7;
options.borderSize = 30;
options.axonScore = 0.3;
options.sourceSize = 2000;
options.recursionSteps = 10;
options.minSize = 100;
options.bboxDist = 1000;
options.voxelSize = [11.24 11.24 28];
options.myelinScore = 0.5;

graph=load('/gaba/u/mberning/results/pipeline/20170217_ROI/graphNewNew.mat')
borderMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalBorder.mat', 'borderSize', 'borderCoM');
segmentMeta = load('/gaba/u/mberning/results/pipeline/20170217_ROI/segmentMeta.mat', 'voxelCount', 'point', 'maxSegId', 'cubeIdx', 'centroid', 'box');
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat', 'p');
segmentMeta = connectEM.addSegmentClassInformation(p, segmentMeta);
globalSegmentPCA = load('/gaba/u/mberning/results/pipeline/20170217_ROI/globalSegmentPCA.mat', 'covMat');
[graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
graph.neighProb = cellfun(@(x)graph.prob(x), neighboursIdx, 'uni', 0);
graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);

segDendrites = arrayfun(@(x)x.nodes(:,4),dendrites,'uni',0);

% Indices of whole cells
idx = [30    16    44    92    53    15    36     4    63    35    33    31    23     8    11    19    13    41    21    29    22    24    12     3    10    7];
versionMapping = [31182 21658 4422 62371 5962 31948 6080 31948 51596 10507 67731 63975 56872 23548 71793 32339 54895 35415 45307 75813 13553 2472 127899 23912 85261 9217]';
idx = versionMapping;

wholeCellDendrites = dendrites(idx);
wholeCellDendriteSegments = segDendrites(idx);

directionality = connectEM.calculateDirectionalityOfAgglomerates(wholeCellDendriteSegments, graph,...
    segmentMeta, borderMeta, globalSegmentPCA, options);
    
% save('/u/mberning/results/pipeline/20170217_ROI/aggloState/dendriteEndings/2017-10-13_directionalityDendrites_1um.mat','directionality')
save('/u/mberning/results/pipeline/20170217_ROI/aggloState/dendriteEndings/2017-10-13_directionalityDendrites_v03_v2_um.mat','directionality')


%%
options.latentScore = 0.5;
options.segDirScore = 0.9;
options.border = [3000; -3000];
options.bboxDist = 1000;

segDendrites = arrayfun(@(x)x.nodes(:,4),dendrites,'uni',0);
segDendrites_v2 = arrayfun(@(x)x.nodes(:,4),dendrites_v2,'uni',0);

directionality = connectEM.calculateDirectionalityOfAgglomerates(segDendrites, graph,...
    segmentMeta, borderMeta, globalSegmentPCA, options);

directionality_v2 = connectEM.calculateDirectionalityOfAgglomerates(segDendrites_v2, graph,...
    segmentMeta, borderMeta, globalSegmentPCA, options);


% idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, directionality.latent, 'uni', 0);
idxDirectional = cellfun(@(x)true(size(x,1),1),directionality.latent,'uni',0);
idxEnding = cellfun(@(x)abs(x) > options.segDirScore, directionality.scores, 'uni', 0);
idxAll = cellfun(@(x,y)find(x&y), idxDirectional, idxEnding, 'uni', 0);
    
thisBorderIdx = cellfun(@(x,y)x(y), directionality.borderIdx, idxAll,'uni',0);
borderPositions = cellfun(@(x) borderMeta.borderCoM(x,:), thisBorderIdx,'uni',0);
 
distance_cutoff = 800;
T = cellfun(@(x) AxonEndings.clusterSynapticInterfaces(x,distance_cutoff,[11.24 11.24 28]), borderPositions, 'uni', 0);


idxDirectional = cellfun(@(x)true(size(x,1),1),directionality_v2.latent,'uni',0);
idxEnding = cellfun(@(x)abs(x) > options.segDirScore, directionality_v2.scores, 'uni', 0);
idxAll = cellfun(@(x,y)find(x&y), idxDirectional, idxEnding, 'uni', 0);
    
thisBorderIdx = cellfun(@(x,y)x(y), directionality_v2.borderIdx, idxAll,'uni',0);
borderPositions = cellfun(@(x) borderMeta.borderCoM(x,:), thisBorderIdx,'uni',0);
 
T_v2 = cellfun(@(x) AxonEndings.clusterSynapticInterfaces(x,distance_cutoff,[11.24 11.24 28]), borderPositions, 'uni', 0);

sum(cell2mat(cellfun(@(x)max(x),T,'uni',0)))
sum(cell2mat(cellfun(@(x)max(x),T_v2,'uni',0)))
