function findEdgesAndBordersWrapper(segFile, edgeFile, borderFile, segmentFile)
    % Computation of edges and borders optimized for 512x512x256

    % Load *SMALL* cube of global segmentation IDs
    load(segFile);

    % Use new function from SynEM for calculation of svg data
    [edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg);

    % Save to files: currently segmentation overwrites old one (now leaves are merged)
    Util.save(edgeFile, edges);
    Util.save(borderFile, borders);
    Util.save(segmentFile, segments);
end
