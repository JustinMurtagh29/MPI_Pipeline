% Author: 
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
% connectEM is run on the graph
% Hierarchical agglomeration is run on the graph
% Both are run here conservatively and the agglomerations from both 
% methods are here combined to get a larger, error-free agglomeration
% state

info = Util.runInfo();

sizeThreshold = 1e6;

% Load parameter
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
load(fullfile(rootDir,'allParameter.mat'))
maxSegId = Seg.Global.getMaxSegId(p);
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
graph = load([p.saveFolder 'graph.mat']);
edges = graph.edges;
Util.log('Loading NC:')
m = load(fullfile(rootDir,'20190405T151728_agglomeration_NCHC', 'NC','NC_20190405T151728_agglos.mat'));
agglosNC = m.agglos;
agglosNCSizes = m.agglosSize;
agglosNC = agglosNC(2:end);

% Keep only agglomerates that have at least sizeThreshold million voxels
idx = agglosNCSizes > sizeThreshold;
agglosNC = agglosNC(idx);
agglosNCSizes = agglosNCSizes(idx); 
clear idx
Util.log('Loading HC:')
m = load(fullfile(rootDir,'20190405T151728_agglomeration_NCHC', 'HC','HC_20190405T151728_agglos.mat'));
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

% do the merging
cc = Graph.findConnectedComponents(edgesBetweenAgglos);
agglosNew = cellfun(@(x)cat(1,agglosHC{x}), cc, 'uni', 0); %#1

% add corresponding aggloNC to new agglos
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
agglosNewAppened = cellfun(@(x,y) unique(cat(1, x, cat(1, outClass{y}))), agglosNew, toAdd, 'uni',0); %#2

% add single HC found by NC but not contained in CC
idxOut = cellfun(@(x) numel(x)==1 & any(x~=0), idxClass);
outClass= idxClass(idxOut);
outClass= cell2mat(outClass)';
outAgglos = agglosNC(idxOut);
notFound = ~ismember(outClass, cat(1,cc{:}));
outClass = outClass(notFound);
outAgglos = outAgglos(notFound);

agglosToAdd = arrayfun(@(x, y) unique(cat(1, outAgglos{x}, agglosHC{y})), [1:numel(outAgglos)]', outClass, 'uni', 0); %#3

% remaining NC
idxOut = cellfun(@(x) numel(x) ==1 & all(x==0), idxClass);
moreAgglos = agglosNC(idxOut); %#4

% remaining HC
idx = setdiff(1:size(agglosHC,1), cat(1,idxClass{:}));
moreMoreAgglos = agglosHC(idx); %#5

% merging all to one
agglosFinal = cat(1, agglosNewAppened, agglosToAdd, moreAgglos, moreMoreAgglos);
agglos = agglosFinal;
agglosSize = cellfun(@(x) sum(segmentMeta.voxelCount(x)), agglos);
[agglosSize, idx] = sort(agglosSize, 'descend');
agglos = agglos(idx);
clear idx

outputFolder = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/20190405T151728_agglomeration_NCHC/simpleMerge/';
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
tic;
Superagglos.skeletonFromAgglo(edges, segmentMeta, ...
    agglosOut, 'agglos', outputFolder, parameters);
toc;















