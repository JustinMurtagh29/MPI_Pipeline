function [graphR, comR, segIdsR, segIdsB] = restrictSGtoRegion(p, graph, com, bbox, excludeIdx)
% Restrict graph and CoMs to bounding box, after removing vessel & somata
% return segIds in bbox region and at border of region as well

    % Load segmentation
    seg = readKnossosRoi(p.seg.root, p.seg.prefix, bbox, 'uint32', '', 'raw'); 
    % Get unique segmentation ids in ROI
    segIdsR = unique(seg(:));
    segIdsR = segIdsR(segIdsR ~= 0);
    segIdsR = segIdsR(~ismember(segIdsR, excludeIdx));
    % Get unique segmentation ids at border
    segBorder = cat(1,reshape(seg([1 end],:,:), [2*size(seg,2)*size(seg,3) 1 1]),...
                      reshape(seg(:,[1 end],:), [2*size(seg,1)*size(seg,3) 1 1]),...
                      reshape(seg(:,:,[1 end]), [2*size(seg,1)*size(seg,2) 1 1]));
    segIdsB = unique(segBorder);
    segIdsB = segIdsB(segIdsB ~= 0);
    segIdsB = segIdsB(~ismember(segIdsB, excludeIdx));
    % Restrict CoMs to bounding box
    idx = false(size(com,1),1);
    idx(segIdsR) = 1;
    comR = com(idx,:);
    % Restrict edges to bounding box
    edgeIdx = ismember(graph.edges(:,1), segIdsR) & ismember(graph.edges(:,2), segIdsR);
    graphR.edges = graph.edges(edgeIdx,:);
    graphR.prob = graph.prob(edgeIdx);
    graphR.cubeLI = graph.cubeLI(edgeIdx);
    graphR.borderCentroid = graph.borderCentroid(edgeIdx,:);

end

