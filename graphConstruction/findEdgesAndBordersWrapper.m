function findEdgesAndBordersWrapper(segFile, edgeFile, borderFile, segmentFile)
    % Computation of edges and borders optimized for 512x512x256

    % Load *SMALL* cube of global segmentation IDs
    load(segFile);

    % Use new function from SynEM for calculation of svg data
    [edges, borders, segments] = SynEM.Svg.findEdgesAndBorders(seg);

    % Save to files: currently segmentation overwrites old one (now leaves are merged)
    save(edgeFile, 'edges');
    save(borderFile, 'borders');
    save(segmentFile, 'segments', '-v7.3');
end
