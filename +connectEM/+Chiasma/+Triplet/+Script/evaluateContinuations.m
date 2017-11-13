% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
nmlDir = '/mnt/mpibr/data/Personal/mottaa/L4/2017-11-09-Finding-Triplets/continuation-tracings';
outDir = '/home/amotta/Desktop/continuation-overlaps';

% based on hard-coded value from `connectEM.generateAxonQueries`
axonFile = fullfile(rootDir, 'aggloState', 'axons_04.mat');

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

axons = load(axonFile, 'axons', 'indBigAxons');
axons = axons.axons(axons.indBigAxons);

%% building look-up table
axonAgglos = Superagglos.getSegIds(axons);
axonLUT = Agglo.buildLUT(maxSegId, axonAgglos);
axonLUT = cat(1, 0, reshape(axonLUT, [], 1));

%% find NML files
nmlFiles = NML.findFiles(nmlDir);
nmlFiles = fullfile(nmlDir, nmlFiles(:));
nmlFileCount = numel(nmlFiles);

%%
nmlNodes = cell(nmlFileCount, 1);
for curIdx = 1:nmlFileCount
    curNmlFile = nmlFiles{curIdx};
    curNml = slurpNml(curNmlFile);
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes.coord = curNodes.coord + 1;
    
    curNodes.segId = Seg.Global.getSegIds(param, curNodes.coord);
    curNodes.aggloId = axonLUT(1 + curNodes.segId);
    
    nmlNodes{curIdx} = curNodes;
end

%% export NML with overlaps
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

for curIdx = 1:nmlFileCount
    curSkel = skeleton();
    
    curNodes = nmlNodes{curIdx};
    curSkel = curSkel.addTree('Manual continuation', curNodes.coord);
    
    curAggloIds = setdiff(curNodes.aggloId, 0);
    for curAggloIdx = 1:numel(curAggloIds)
        curAggloId = curAggloIds(curAggloIdx);
        curAgglo = axons(curAggloId);
        
        curSkel = curSkel.addTree( ...
            sprintf('Agglomerate %d', curAggloId), ...
            curAgglo.nodes(:, 1:3), curAgglo.edges);
    end
    
    
    curFileName = nmlFiles{curIdx};
   [~, curFileName] = fileparts(curFileName);
    curFileName = fullfile(outDir, strcat(curFileName, '.nml'));
   
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel.write(curFileName);
end