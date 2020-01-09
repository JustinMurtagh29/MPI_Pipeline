function skel = inspectTPs(param, curCount, className, gt)

curPoints = Seg.Global.getSegToPointMap(param);

curMask = gt.class == className;

tp = table;
tp.segId = gt.segId;
tp.label = gt.label(:, curMask);
tp.score = gt.scores(:, curMask);
tp.probs = gt.probs(:,curMask);

tp = tp(tp.label < 0, :);
tp = sortrows(tp, 'probs', 'descend');

curCount = min(curCount,sum(tp.label < 0));
curDigits = ceil(log10(1 + curCount));

skel = skeleton();
skel = skel.setParams(param.experimentName,param.raw.voxelSize,[1 1 1]);
curSegIds = tp.segId(1:curCount);
curScores = tp.score(1:curCount);
curProbs = tp.probs(1:curCount);

curNodes = curPoints(curSegIds, :);
curNames = arrayfun( ...
    @(idx, segId, prob, score) sprintf( ...
        '%0*d. Segment %d. Prob %.2f. Score %.3f', ...
        curDigits, idx, segId, prob, score), ...
    (1:curCount)', curSegIds, curProbs, curScores, ...
    'UniformOutput', false);

skel = skel.addNodesAsTrees(curNodes, curNames);
end
