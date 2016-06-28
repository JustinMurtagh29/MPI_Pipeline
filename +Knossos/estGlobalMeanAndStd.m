function [meanVal, stdVal] = estGlobalMeanAndStd(p)
    % [meanVal, stdVal] = estGlobalMeanAndStd(p)
    %   Estimates the global mean and standard deviation of
    %   the raw data by randomly sampling from small cubes.
    %
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>

    % config
    nrCubesToSample = 50;
    cubeSize = 100;

    % How many 100^3 samples to use for determination of normalization
    % constants, 100 seems rather too much but ok as only takes 5 min
    sizeOfRoi = p.bbox(:, 2) - p.bbox(:, 1) + 1;
    meanVal = zeros(nrCubesToSample, 1);
    stdVal = zeros(nrCubesToSample, 1);

    for i=1:nrCubesToSample
        % randomly pick lower left corner
        lowerLeft = [ ...
            randi(sizeOfRoi(1) - (cubeSize - 1)); ...
            randi(sizeOfRoi(2) - (cubeSize - 1)); ...
            randi(sizeOfRoi(3) - (cubeSize - 1))];
        lowerLeft = lowerLeft + p.bbox(:, 1) - 1;

        % build bounding box
        bbox = [ ...
            lowerLeft, ...
            lowerLeft + (cubeSize - 1)];

        % load raw data
        raw = loadRawData( ...
            p.raw.root, p.raw.prefix, bbox);

        % compute local mean and std
        [meanVal(i), stdVal(i)] = ...
            Knossos.calcMeanAndStd(raw);
    end
    
    if any(meanVal == 0) || any(stdVal == 0)
        nZero = max(length(find(meanVal == 0)), length(find(stdVal == 0)));
        warning('Found %d of %d randomly chosen cubes with 0 mean or standard deviation in bounding box', nZero, nrCubesToSample);
    end

    meanVal = median(meanVal);
    stdVal = median(stdVal);
end
