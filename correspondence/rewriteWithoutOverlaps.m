function rewriteWithoutOverlaps( inFile, outFile, bbS, bbB)

    segOld = extractBboxSmall(inFile, bbS, bbB);
    % Recalculate connected components of segments in bounding box small
    seg = uint32(bwlabeln(segOld > 0, 26));
    [oldSegments, newSegments] = determineChanges(segOld, seg);
    numEl = calcNumberSegments(seg);
    Util.save(outFile, seg, numEl, oldSegments, newSegments);

end

function [segIdsOld, segIdsNew] = determineChanges(old, new)
    uniqueRows = unique(cat(2, uint32(old(:)), new(:)), 'rows');
    uniqueRows = uniqueRows(~any(uniqueRows == 0, 2),:);
    segIdsOld = uniqueRows(:,1);
    segIdsNew = uniqueRows(:,2);
end

