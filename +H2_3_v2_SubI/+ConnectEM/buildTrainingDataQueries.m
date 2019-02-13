% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

boxesWk = [ ...
    9120, 6595, 2982, 267, 267, 108; ...
    11599, 7708, 2779, 267, 267, 108];
sampleRange = [402, 600];

info = Util.runInfo();
Util.showRunInfo(info);
%% Configuration
rootDir = '/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet_HC/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

outDir = fullfile(param.saveFolder,'tracings','connectEM');

segPoints = Seg.Global.getSegToPointMap(param);

% Build graph
graph = table;

curData = fullfile(param.saveFolder, 'globalEdges.mat');
curData = load(curData, 'edges');
graph.edge = curData.edges;
clear curData;

curData = fullfile(param.saveFolder, 'globalBorder.mat');
curData = load(curData, 'borderCoM');
graph.com = curData.borderCoM;
clear curData;

%% Sample edges
clear cur*:

for curBoxIdx = 1:size(boxesWk, 1)
    curBoxWk = boxesWk(curBoxIdx, :);
    curBox = Util.convertWebknossosToMatlabBbox(curBoxWk);
    
    curSegIds = loadSegDataGlobal(param.seg, curBox);
    curSegIds = reshape(setdiff(curSegIds, 0), [], 1);
    
    curBorderIds = ismember(graph.edge, curSegIds);
    curBorderIds = find(all(curBorderIds, 2));
    
    curRange = min(sampleRange, numel(curBorderIds));
    curRange = curRange(1):curRange(2);
    
    rng(0);
    curRandIds = randperm(numel(curBorderIds));
    curRandIds = curBorderIds(curRandIds(curRange));
    
    curDigits = ceil(log10(1 + numel(curRandIds)));
    
    curSkel = skeleton();
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);
    
    for curBorderIdx = 1:numel(curRandIds)
        curBorderId = curRandIds(curBorderIdx);
        
        curName = sprintf( ...
            '%0*d. Border %d', curDigits, ...
            curBorderIdx, curBorderId);
        
        curNodes = [ ...
            double(segPoints(graph.edge(curBorderId, 1), :)); ...
            double(graph.com(curBorderId, :)); ...
            double(segPoints(graph.edge(curBorderId, 2), :))];
        
        curSkel = curSkel.addTree(curName, curNodes);
    end
    
    if ~isempty(outDir)
        curSkel.write(fullfile(outDir, sprintf('box-%d.nml', curBoxIdx)));
    end
end
