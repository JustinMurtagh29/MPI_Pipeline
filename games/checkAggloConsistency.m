function checkAggloConsistency(agglo)
% Check whether supervoxel graph agglomeration is consitent
% Zero ID and Overlap in Components currenlty checked

    labels = fieldnames(agglo);
    for i=1:length(labels)
        segIds{i} = cat(1,agglo.(labels{i}){:});
        % Check for background in all components
        zeroHit = any(segIds{i} == 0);
        if zeroHit
            display(['Zero segment ID found in ' labels{i}]);
        end
        twicePresent = length(segIds{i}) ~= length(unique(segIds{i}));
        if twicePresent
            display(['Segment ID present twice in ' labels{i}]);
            nrCells = length(agglo.(labels{i}));
            for j=1:nrCells
                for k=(j+1):nrCells
                    ids = intersect(agglo.(labels{i}){j},agglo.(labels{i}){k});
                    if length(ids) > 0
                        display(['Segment ID(s) present in ' num2str(j) ' and ' num2str(k) ':']);
                        display(num2str(ids));
                    end
                end
            end
        end
    end
    
    % Check for duplicated between labels
    for i=1:length(segIds)
        for j=(i+1):length(segIds)
            overlapPresent = length(intersect(segIds{i},segIds{j})) > 0;;
            if overlapPresent
                display(['Segment ID overlap between ' labels{i} ' and ' labels{j}]);
            end
        end
    end

end

