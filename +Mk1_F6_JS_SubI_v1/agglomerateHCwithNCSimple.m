% Author: 
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
% connectEM is run on the graph
% Hierarchical agglomeration is run on the graph
% Both are run here conservatively and the agglomerations from both 
% methods are here combined to get a larger, error-free agglomeration
% state

info = Util.runInfo();

sizeThreshold = 1e6;
datasetName= 'Mk1_F6_JS_SubI_v1';
timeStamp = '20190409T105848';
% Load parameter
rootDir = ['/tmpscratch/sahilloo/data/' datasetName '/pipelineRun_mr2e_wsmrnet/'];
load(fullfile(rootDir,'allParameter.mat'))
maxSegId = Seg.Global.getMaxSegId(p);
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
graph = load([p.saveFolder 'graph.mat']);
graphHC = load(fullfile(rootDir, 'agglomeration/graph.mat'));
graphHC = graphHC.graph;
edges = unique(cat(1,graphHC.edges, graph.edges),'rows','stable');

Util.log('Loading NC:')
m = load(fullfile(rootDir,[timeStamp, '_agglomeration_NCHC'], 'NC',['NC_' timeStamp '_agglos.mat']));
agglosNC = m.agglos;
agglosNCSizes = m.agglosSize;
agglosNC = agglosNC(2:end);

% Keep only agglomerates that have at least sizeThreshold million voxels
idx = agglosNCSizes > sizeThreshold;
agglosNC = agglosNC(idx);
agglosNCSizes = agglosNCSizes(idx); 
clear idx
Util.log('Loading HC:')
m = load(fullfile(rootDir,[timeStamp '_agglomeration_NCHC'], 'HC',['HC_' timeStamp '_agglos.mat']));
agglosHC = m.agglos;
agglosHC = agglosHC(2:end);

Util.log('Combining without individual mega mergers:')
% do the thing
lut = Agglo.buildLUT(maxSegId, agglosHC);
edgesBetweenAgglos = [];
idxClass = {};
for i=1:size(agglosNC,1)
    curAgglo = agglosNC{i};
    idx = unique(nonzeros(lut(curAgglo)));
    if numel(idx)>1
        curEdges = cat(2,idx(1:end-1), idx(2:end));
        edgesBetweenAgglos = cat(1, edgesBetweenAgglos, curEdges);
    end
    idxClass{i} = idx;
end

Util.log('do the merging')
tic;
cc = Graph.findConnectedComponents(edgesBetweenAgglos); % keeps only agglos above size 1
agglosNew = cellfun(@(x) cat(1,agglosHC{x}), cc, 'uni', 0); %#1
toc;

Util.log('add corresponding aggloNC that contributed to cc above to new agglos')
tic;
idxOut = cellfun(@(x) numel(x)>1, idxClass);
agglosOut = agglosNC(idxOut);
outClass = idxClass(idxOut)';
maxId = max(cat(1,outClass{:}));
lutOut = Agglo.buildLUT(maxId, outClass);
toAdd ={};
for i=1:size(cc,1)
    curCC = cc{i};
    toAdd{i} = unique(nonzeros(lutOut(curCC)));
end
toAdd = toAdd';
agglosNewAppended = cellfun(@(x,y) unique(cat(1, x, cat(1, agglosOut{y}))),agglosNew, toAdd, 'uni',0); %#2
toc;

Util.log('add corresponding aggloNC that have been found by other cc to new appended agglos')
tic;
idxOut = cellfun(@(x) numel(x)==1 & any(x~=0), idxClass);
agglosOut = agglosNC(idxOut);
outClass = idxClass(idxOut);
outClass = cell2mat(outClass);
clear toAdd
toAdd = cellfun(@(x) ismember(outClass',x),cc, 'uni',0);
agglosNewAppended2 = cellfun(@(x,y) unique(cat(1, x, cat(1, agglosOut{y}))), agglosNewAppended, toAdd, 'uni',0); %#3
toc;

Util.log('add single HC found by NC but not contained in CC')
notFound = ~ismember(outClass, cat(1,cc{:}));
outClassNotFound = outClass(notFound);
agglosOutNotFound = agglosOut(notFound);
outClassNotFoundUnique = unique(outClassNotFound);
agglosToAdd1 = arrayfun(@(x) unique(cat(1, agglosHC{x}, agglosOutNotFound{outClassNotFound==x} )),outClassNotFoundUnique, 'uni', 0); %#3

% remaining NC
idxOut = cellfun(@(x) numel(x) ==1 & all(x==0), idxClass);
agglosToAdd2 = agglosNC(idxOut); %#4

% remaining HC
idx = setdiff(1:size(agglosHC,1), cat(1,idxClass{:}));
agglosToAdd3 = agglosHC(idx); %#5

Util.log('merging all to one')
agglosFinal = cat(1, agglosNewAppended2(:), agglosToAdd1(:), agglosToAdd2(:), agglosToAdd3(:));
agglos = agglosFinal;
agglosSize = cellfun(@(x) sum(segmentMeta.voxelCount(x)), agglos);
[agglosSize, idx] = sort(agglosSize, 'descend');
agglos = agglos(idx);
clear idx

outputFolder = fullfile(rootDir, [timeStamp '_agglomeration_NCHC'], 'simpleMerge');
mkdir(outputFolder)
agglosOut = agglos(1:100);
display('Writing skeletons for debugging the process:');
parameters.experiment.name= p.experimentName;
parameters.scale.x = num2str(p.raw.voxelSize(1));
parameters.scale.y = num2str(p.raw.voxelSize(2));
parameters.scale.z = num2str(p.raw.voxelSize(3));
parameters.offset.x = '0';
parameters.offset.y = '0';
parameters.offset.z = '0';

edges = unique(cat(1,graphHC.edges, graph.edges),'rows','stable');
tic;
Superagglos.skeletonFromAgglo(edges, segmentMeta, ...
    agglosOut, 'agglos', outputFolder, parameters);
toc;















