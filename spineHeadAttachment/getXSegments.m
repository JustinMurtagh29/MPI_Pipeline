function segments = getXSegments(graph, positiveList, threshold, maxNeighbours, documentID, maxSteps, toexclude)
result = agglomerate_fast(graph, threshold, {}, toexclude, maxNeighbours);
if ~isnan(documentID)
    biggest = find(result.components == result.components(documentID));
    save('documentID', biggest);
end
segments = false(size(result.connM,1), 1);
segments(positiveList(positiveList <= length(segments))) = true;
connM_sym = result.connM + result.connM';
connM_sym(connM_sym > 1) = 1;
counter = 0;
while true
    'step'
    counter = counter + 1;
    segmentsNew = logical(connM_sym * segments);
    if isequal(segmentsNew, segments) || counter > maxSteps
        break
    else
        segments = segmentsNew;
    end
end