function skel = inspectTPs(param, curCount, className, gt)

curPoints = Seg.Global.getSegToPointMap(param);

curMask = gt.class == className;

fp = table;
fp.segId = gt.segId;
fp.label = gt.label(:, curMask);
fp.score = gt.scores(:, curMask);
fp.probs = gt.probs(:,curMask);

fp = fp(fp.label > 0, :);
fp = sortrows(fp, 'probs', 'ascend');

curCount = min(curCount,sum(fp.label>0));
curDigits = ceil(log10(1 + curCount));

skel = skeleton();
skel = skel.setParams(param.experimentName,param.raw.voxelSize,[1 1 1]);
curSegIds = fp.segId(1:curCount);
curScores = fp.score(1:curCount);
curProbs = fp.probs(1:curCount);

curNodes = curPoints(curSegIds, :);
curNames = arrayfun( ...
    @(idx, segId, prob, score) sprintf( ...
        '%0*d. Segment %d. Prob %.2f. Score %.3f', ...
        curDigits, idx, segId, prob, score), ...
    (1:curCount)', curSegIds, curProbs, curScores, ...
    'UniformOutput', false);

skel = skel.addNodesAsTrees(curNodes, curNames);
end
