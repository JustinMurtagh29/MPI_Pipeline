function segIds = determineSegIdsAtSGborders(p, excludeIds)
    % Extract segmentation IDs at each face of supervoxel cube    

    bbox = p.bbox;
    bbox(1,1) = bbox(1,2);
    segIds = unique(readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'));
    bbox = p.bbox;
    bbox(1,2) = bbox(1,1);
    segIds = [segIds; unique(readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'))];
    bbox = p.bbox;
    bbox(2,1) = bbox(2,2);
    segIds = [segIds; unique(readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'))];
    bbox = p.bbox;
    bbox(2,2) = bbox(2,1);
    segIds = [segIds; unique(readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'))];
    bbox = p.bbox;
    bbox(3,1) = bbox(3,2);
    segIds = [segIds; unique(readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'))];
    bbox = p.bbox;
    bbox(3,2) = bbox(3,1);
    segIds = [segIds; unique(readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'))];

    % Remove zeroes
    segIds = unique(segIds(segIds ~= 0));
    % Remove excludeIds (nuclei and vessel IDs in first try)
    segIds(ismember(segIds, excludeIds)) = [];
    % Cast to cell for agglomeration
    segIds = mat2cell(segIds, ones(length(segIds),1), 1);
 
end
