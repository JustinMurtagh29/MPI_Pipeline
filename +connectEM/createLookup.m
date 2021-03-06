function agglos_reverse = createLookup(segmentMeta, agglos)
    % Build lookup of agglomerate ID based on segment ID
    agglos_reverse = zeros(segmentMeta.maxSegId, 1);
    for idx = 1 : length(agglos)
        agglos_reverse(agglos{idx}) = idx;
    end
end

