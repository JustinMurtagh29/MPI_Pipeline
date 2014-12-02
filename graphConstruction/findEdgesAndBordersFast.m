function findEdgesAndBordersFast(seg, edgeFile, borderFile, bbox)
    % Computation of edges and borders optimized for 512x512x256

    % Load segmentation from file
    %[seg, segSmall] = loadSegData(segFile, tileBorder);
    % Load global segmentation instead of local
    segSmall = loadSegDataGlobal(seg.root, seg.prefix, bbox);
    segSmall = int32(segSmall);

    % Pad array with a 1 voxel surround with a new unique value
    globalBorderId = max(segSmall(:))+1;
    segSmall = padarray(segSmall,[1 1 1],globalBorderId);

    % Size of padded segmentation
    [M,N,P] = size(segSmall);

    % Construct 26-connectivity linear indices shift for padded segmentation
    vec = int32([(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])]);

    % Find linear inidices of all wall voxel
    ind = find(segSmall==0);

    % Find segmentation ID of all neighbours of all wall voxel (according to 26
    % connectivity)
    nInd = bsxfun(@plus, ind, vec);
    nSegId = segSmall(nInd);

    % Find edges
    edges = findEdges(nSegId);

    % Remove edges to padded value (see first comment)
    % This uses fact that globalBorderId is largest as is made sure in the beginning 
    edgesToBorder = edges(edges(:,2) == globalBorderId,:);
    edges(edges(:,2) == globalBorderId,:) = [];

    % Find borders
    [edges,borders] = findBorders(ind, nInd, nSegId, edges, M, N, P);

    % Save to files: currently segmentation overwrites old one (now leaves are merged)
    save(edgeFile, 'edges', 'edgesToBorder');
    save(borderFile, 'borders');

end

