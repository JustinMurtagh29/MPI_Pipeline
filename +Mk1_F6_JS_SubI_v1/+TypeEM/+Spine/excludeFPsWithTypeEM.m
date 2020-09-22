% run this after debugBoxOverAndUnder.m

segmentMeta = load([param.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';

segmentPredictions = load([param.saveFolder 'segmentAggloPredictions.mat'], '-mat');
% ... glia
segmentMeta.gliaProb = zeros(segmentMeta.maxSegId, 1);
idx = ~isnan(segmentPredictions.probs(:,1));
segmentMeta.gliaProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,1);
% ... axon
segmentMeta.axonProb = zeros(segmentMeta.maxSegId, 1);
idx = ~isnan(segmentPredictions.probs(:,2));
segmentMeta.axonProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,2);
% ... dendrite
segmentMeta.dendriteProb = zeros(segmentMeta.maxSegId, 1);
idx = ~isnan(segmentPredictions.probs(:,3));
segmentMeta.dendriteProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,3);

% Test set ground-truth true labels in predictions
probThr = 0.50;
idxPred = curGtTest.probs>probThr;
idxMissed = curGtTest.label == 1 & ~idxPred;

predSegIds = curGtTest.segId(idxPred);
sprintf('Predicted SHs above %.2f prob thr: %d', probThr, sum(idxPred))
sprintf('Predicted SHs missed above %.2f prob thr: %d', probThr, sum(idxMissed))

idxFPs = curGtTest.label(idxPred) == -1;
sprintf('FPs in predictions: %d', sum(idxFPs))

idxTPs = ~idxFPs;

minAxonProb = 0.6; % GA on agglomeration for dendrites
idxAxonFPs = segmentMeta.axonProb(predSegIds(idxFPs)) > minAxonProb;
sprintf('FPs excluded with minAxonProb %.2f : %d', minAxonProb, sum(idxAxonFPs))

idxAxonTPs = segmentMeta.axonProb(predSegIds(idxTPs)) > minAxonProb;
sprintf('TPs excluded with minAxonProb %.2f : %d', minAxonProb, sum(idxAxonTPs))



