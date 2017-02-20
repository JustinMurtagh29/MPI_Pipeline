function [segmentOverlap, uniqueSegId] = getSegmentHeuritsicsScore(seg, binaryMap, bbox)
    % Takes segmentation struct (root, prefix, bbox) and looks up the fraction 
    % of each segment overlapping with binaryMap (specified by root, prefix and segId to lookup)

    seg = loadSegDataGlobal(seg, bbox);
    uniqueSegId = unique(seg(seg ~= 0));
    for i=1:length(binaryMap)
        class = loadSegDataGlobal(binaryMap(i), bbox);
        class = class == binaryMap(i).segId;
        segOccurences = histc(seg(seg~=0), uniqueSegId); 
        classVoxels = histc(seg(seg~=0 & class), uniqueSegId);
        segmentOverlap{i} = classVoxels ./ segOccurences;       
    end

end

