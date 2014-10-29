function trueEdges = approximateGroundTruth(seg, seg_manuell)

% Find correspondence of manuell and automated segmentation
a = unique(seg_manuell);
a(a == 0) = [];
for i=1:length(a)
    withinSegment{i} = seg(seg_manuell == a(i));
    withinUnique{i} = unique(withinSegment{i});
    withinUnique{i}(withinUnique{i} == 0) = [];  
    counts{i} = histc(withinSegment{i},withinUnique{i}); 
    withinUnique{i}(counts{i} < 100) = [];
end

nrIDs = length(unique(seg));
nrIDs(nrIDs == 0) = [];
trueEdges = [];
for i=1:size(withinUnique,2)
    for j=1:size(withinUnique{i},1)
        for k=1:size(withinUnique{i},1)
            if k~=j
                trueEdges(end+1,:) = [withinUnique{i}(k) withinUnique{i}(j)];
            end
        end
    end
end
trueEdges = sort(trueEdges,2);
trueEdges = unique(trueEdges, 'rows');

end

