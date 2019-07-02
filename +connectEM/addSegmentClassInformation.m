function segmentMeta = addSegmentClassInformation(p, segmentMeta)
    % This function loads the segment class probabilties stored by Alessandro
    % It creates seperate fields in the structure for each segment class and changes to indexing vector
    % Convenience funtion

    % Use agglomerate based statistics for everything else
    segmentPredictions = load([p.saveFolder 'segmentAggloPredictions.mat'], '-mat');
    % ... glia
    segmentMeta.gliaProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probsMulti(:,1));
    segmentMeta.gliaProb(segmentPredictions.segId(idx)) = segmentPredictions.probsMulti(idx,1);
    % ... axon
    segmentMeta.axonProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probsMulti(:,2));
    segmentMeta.axonProb(segmentPredictions.segId(idx)) = segmentPredictions.probsMulti(idx,2);
    % ... dendrite
    segmentMeta.dendriteProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probsMulti(:,3));
    segmentMeta.dendriteProb(segmentPredictions.segId(idx)) = segmentPredictions.probsMulti(idx,3);
    % ... spine head (not multi-class)
%{
    segmentPredictions = load([p.saveFolder 'segmentPredictions.mat']);
    segmentMeta.spineProb = zeros(segmentMeta.maxSegId, 1);
    idx = ~isnan(segmentPredictions.probs(:,4));
    segmentMeta.spineProb(segmentPredictions.segId(idx)) = segmentPredictions.probs(idx,4);
%}
end

