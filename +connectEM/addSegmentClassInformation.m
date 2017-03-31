function segmentMeta = addSegmentClassInformation(p, segmentMeta)
    % This function loads the segment class probabilties stored by Alessandro
    % It creates seperate fields in the structure for each segment class and changes to indexing vector
    % Convenience funtion

    segmentPredictions = load([p.saveFolder 'segmentPredictions.mat']);
    segmentPredictions.probs(:,1:3) = bsxfun(@rdivide, segmentPredicitions.probs(:,1:3), sum(segmentPredictions.probs(:,1:3),2));
    % ... glia
    segmentMeta.gliaProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probs(:,1));
    segmentMeta.gliaProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,1);
    segmentMeta.isGlia = segmentMeta.gliaProb > 0.5;
    % ... axon
    segmentMeta.axonProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probs(:,2));
    segmentMeta.axonProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,2);
    segmentMeta.isAxon = segmentMeta.axonProb > 0.5;
    % ... dendrite
    segmentMeta.dendriteProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probs(:,3));
    segmentMeta.dendriteProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,3);
    segmentMeta.isDendrite = segmentMeta.dendriteProb > 0.5;
    % ... spine head
    segmentMeta.spineProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probs(:,4));
    segmentMeta.spineProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,4);
    segmentMeta.isSpine = segmentMeta.spineProb > 0.5;

end

