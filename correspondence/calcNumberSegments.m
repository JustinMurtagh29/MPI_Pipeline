function numEl = calcNumberSegments(seg)
    % Calculate offset local to global segmentation for each cube

    % How many unique global IDs are needed for this cube
    globalIds = unique(seg(:));
    nrGlobalIdsNeededForThisCube = length(globalIds);
    % If there is a 0, do not count this one 
    if any(globalIds == 0)
        nrGlobalIdsNeededForThisCube = nrGlobalIdsNeededForThisCube - 1;
    end
    numEl = uint32(nrGlobalIdsNeededForThisCube); 

end

