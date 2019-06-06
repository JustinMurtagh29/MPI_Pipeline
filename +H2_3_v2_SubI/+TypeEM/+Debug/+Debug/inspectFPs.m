function skel = inspectFPs(param, curCount, className, gt)

className = 'glia';
clear cur*;
curCount = 50;
curDigits = ceil(log10(1 + curCount));
curPoints = Seg.Global.getSegToPointMap(param);

curMask = curGtTest.class == className;

fp = table;
fp.segId = gt.segId;
fp.label = gt.label(:, curMask);
fp.score = gt.scores(:, curMask);

fp = fp(fp.label < 0, :);
fp = sortrows(fp, 'score', 'descend');

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

curSegIds = fp.segId(1:curCount);
curScores = fp.score(1:curCount);

curNodes = curPoints(curSegIds, :);
curNames = arrayfun( ...
    @(idx, segId, score) sprintf( ...
        '%0*d. Segment %d. Score %.3f', ...
        curDigits, idx, segId, score), ...
    (1:curCount)', curSegIds, curScores, ...
    'UniformOutput', false);

skel = skel.addNodesAsTrees(curNodes, curNames);
end
