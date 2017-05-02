axonsFinalAll = [axonsFinal ...
    num2cell(setdiff(find(segmentMeta.axonProb > 0.5), cell2mat(axonsFinal)))];
todo = [];
for idx = 1 : length(axonsFinalAll)
    for idx2 = 1 : length(axonsFinalAll{idx})
        % detect local surround

        % calculate PCA of local surround (Alessandro)

        % find all outgoing edges of current segment

        % calculate minmax score of those

        todo = [todo; outgoing & memberminmax];
    end
end
