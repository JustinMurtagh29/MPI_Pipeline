function [segmentOverlap, uniqueSegId,voxelCount] = getSegmentHeuristicsScore(seg, labelMap, bbox)
    % Takes segmentation struct (root, prefix, bbox) and looks up the fraction 
    % of each segment overlapping with labelMap (specified by root, prefix/backend and segId to lookup)
    % if segId has to values, the function uses all values between these
    % values (including the left side)
    % if the third output is defined, it also returns the voxel count of
    % each connected component in the labelMap which is associated with the
    % segment that overlaps the most with the connected component.

    seg = loadSegDataGlobal(seg, bbox);
    uniqueSegId = unique(seg(seg ~= 0));
    segmentOverlap = cell(length(labelMap),1);
    if nargout > 2
        voxelCount = cell(length(labelMap),1);
    end
    for i=1:length(labelMap)
        class = loadSegDataGlobal(labelMap(i), bbox);
        if numel(labelMap(i).segId) == 2
            class = class >= labelMap(i).segId(1) & class < labelMap(i).segId(2); % make logical map with all values between segId(1) and segId(2)
        else
            class = class == labelMap(i).segId; % make logical map with all values equal to segId
        end
        segOccurences = histc(seg(seg~=0), uniqueSegId); 
        classVoxels = histc(seg(seg~=0 & class), uniqueSegId);
        if size(classVoxels,2) > 1
            % If only one voxel is found, histc switches dimensionality of output for some reason
            segmentOverlap{i} = classVoxels' ./ segOccurences;
        else
            segmentOverlap{i} = classVoxels ./ segOccurences;
        end
        if nargout > 2
            % calculate the voxelCount of each labeled connected component (CC) and associate it
            % with the segment it overlaps the most
            class = bwlabeln(class); % label each CC with a different color
            voxelCounts = histc(class(:),1:max(class(:))); % calculate the voxelCount of each CC
            mostOverlapSegs = accumarray(class(class~=0 & seg~=0),seg(class~=0 & seg~=0),[],@mode); % calculate which segIds overlap the most with each CC
            % now treat cases were a label only lies on border and thus
            % could not be counted. first be sure that mostOverlapSegs has
            % the same size as the number of labels
            if numel(mostOverlapSegs) < numel(voxelCounts)
                mostOverlapSegs(numel(voxelCounts)) = 0;
            end
            % now remove not found labels (because they were exclusively on the
            % segment borders)
            voxelCounts(mostOverlapSegs==0) = [];
            mostOverlapSegs(mostOverlapSegs==0) = [];
            voxelCount{i} = zeros(numel(uniqueSegId),1);
            [~,ind] = ismember(mostOverlapSegs,uniqueSegId);
            voxelCount{i}(ind) = voxelCounts;  % associate these segIds with the voxelCount
        end
    end

end

