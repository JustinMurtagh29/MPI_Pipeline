function curSkel = inspectFPs(param, gt, scores, scoreThr )

info = Util.runInfo();
Util.showRunInfo(info);

% write out FPs above the scoreThr from the gt
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

clear cur*

% find FPs
idxAbove = scores > scoreThr;
borderIdsAbove = gt.borderId(idxAbove);
labelsAbove = gt.label(idxAbove);

% ignore double borders for inspection
[borderIdsAbove, idxUnique] = unique(borderIdsAbove);
labelsAbove = labelsAbove(idxUnique);

idxFP = ~(labelsAbove == 1);
curBorderIds  = borderIdsAbove(idxFP);

curSkel = skeleton();
curSkel = Skeleton.setParams4Pipeline(curSkel, param);
curSkel = Skeleton.setDescriptionFromRunInfo(curSkel, info);

for curBorderIdx = 1:numel(curBorderIds)
    curBorderId = curBorderIds(curBorderIdx);
    
    curName = sprintf( ...
        'Border %d', curBorderId);
    
    curNodes = [ ...
        double(segPoints(graph.edge(curBorderId, 1), :)); ...
        double(graph.com(curBorderId, :)); ...
        double(segPoints(graph.edge(curBorderId, 2), :))];
    
    curSkel = curSkel.addTree(curName, curNodes);
end

end





