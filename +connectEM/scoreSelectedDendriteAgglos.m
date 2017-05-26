function metrics = scoreSelectedDendriteAgglos(segmentMeta, dendritesFinal, p)

counter = 0;
for idx = 1 : length(dendritesFinal)
    volume = sum(segmentMeta.voxelCount(dendritesFinal{idx}));
    scores = sort(segmentMeta.dendriteProb(dendritesFinal{idx}));
    if volume > 100000 && scores(round(end*3/4)) > 0.5
        counter = counter + 1;

        thissize = max(pdist(bsxfun(@times, segmentMeta.point(:, dendritesFinal{idx})', [11.24, 11.24, 28])));
        metrics(counter).idx = idx;
        metrics(counter).volsize = thissize;
        allboxes = transpose(reshape(segmentMeta.box(:, :, dendritesFinal{idx}), 3, []));
        mybbox = [p.bbox(1,1), Inf, Inf;
                  p.bbox(1,2), Inf, Inf;
                  Inf, p.bbox(2,1), Inf;
                  Inf, p.bbox(2,2), Inf;
                  Inf, Inf, p.bbox(3,1);
                  Inf, Inf, p.bbox(3,2)];
        touching = pdist2(mybbox, allboxes, @(x,y)min(bsxfun(@times, abs(bsxfun(@minus,x,y)), [11.24, 11.24, 28]),[], 2)) < 2000;
        metrics(counter).singletouch = any(touching(:));
        metrics(counter).doubletouch = any(sum(touching) > 1) && any(sum(touching, 2) > 1);
    end
end
