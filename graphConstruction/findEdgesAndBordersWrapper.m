function findEdgesAndBordersWrapper(segFile, edgeFile, borderFile, segmentFile, excludeVoxels)
    % Computation of edges and borders optimized for 512x512x256

    % Load *SMALL* cube of global segmentation IDs
    load(segFile);

    % Use new function from SynEM for calculation of svg data
    if nargin > 4
        if all(excludeVoxels(:)) % skip calculation if completely in mirrorPaded region
            return
        end
        [edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg,excludeVoxels);
    else
        [edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg);
    end

    % Save to files: currently segmentation overwrites old one (now leaves are merged)
    Util.save(edgeFile, edges);
    Util.save(borderFile, borders);
    Util.save(segmentFile, segments);
end
