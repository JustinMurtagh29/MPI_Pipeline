function findEdgesAndBordersFast(segFile, edgeFile, borderFile, segmentFile)
    % Computation of edges and borders optimized for 512x512x256

    % Load *SMALL* cube of global segmentation IDs
    load(segFile);
    segSmall = int32(seg);

    % Pad array with a 1 voxel surround with a new unique value
    globalBorderId = max(segSmall(:))+1;
    seg = padarray(segSmall,[1 1 1],globalBorderId);

    % Size of padded segmentation
    [M,N,P] = size(seg);

    % Construct 26-connectivity linear indices shift for padded segmentation
    vec = int32([(-M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1]) [-M-1 -M -M+1 -1 1 M-1 M M+1] (M*N+[-M-1 -M -M+1 -1 0 1 M-1 M M+1])]);

    % Find linear inidices of all wall voxel
    ind = int32(find(seg==0));

    % Find segmentation ID of all neighbours of all wall voxel (according to 26
    % connectivity)
    nInd = bsxfun(@plus, ind, vec);
    nSegId = seg(nInd);

    % Find edges
    edges = findEdges(nSegId);

    % Remove edges to padded value (see first comment)
    % This uses fact that globalBorderId is largest as is made sure in the beginning 
    edgesToBorder = edges(edges(:,2) == globalBorderId,:);
    edges(edges(:,2) == globalBorderId,:) = [];

    % Find borders
    [edges,borders] = findBorders(ind, nInd, nSegId, edges, M, N, P);

    % sort edges and borders
    [edges, sortedRows] = sortrows(edges);
    borders = borders(sortedRows);

    % Save to files: currently segmentation overwrites old one (now leaves are merged)
    save(edgeFile, 'edges', 'edgesToBorder');
    save(borderFile, 'borders');

    % Get segment pixelIdxLists
    segments = regionprops( ...
        segSmall, segSmall, 'PixelIdxList', 'MinIntensity');
    
    % NOTE
    %   Previously, regionprops was run on the local segmentation volume
    %   where the segment IDs started with one. Now, we have to skip
    %   all IDs below the smallest occuring global segment ID.
    segments = segments(arrayfun( ...
        @(curSeg) ~isempty(curSeg.MinIntensity), segments));
    
    % NOTE
    %   We use the mean intensity to find the global segment ID.
    [segments.Id] = segments.MinIntensity;
    
    % remove MinIntensity field
    segments = rmfield(segments, 'MinIntensity');
    save(segmentFile, 'segments', '-v7.3');
end

