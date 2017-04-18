function [y, ysize] = scoreSelectedDendriteAgglos(segmentMeta, dendritesFinal)
y = [];
ysize = [];
for idx = 1 : length(dendritesFinal)
    volume = sum(segmentMeta.voxelCount(dendritesFinal{idx}));
    scores = sort(segmentMeta.dendriteProb(dendritesFinal{idx}));
    if volume > 100000 && scores(round(end*3/4)) > 0.5
        thissize = max(pdist(bsxfun(@times, segmentMeta.point(:, dendritesFinal{idx})', [11.24, 11.24, 28])));
        y = [y, idx];
        ysize = [ysize, thissize];
    end
end
