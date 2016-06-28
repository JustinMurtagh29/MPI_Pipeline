function [meanVal, stdVal] = calcMeanAndStd(raw)
    % [meanVal, stdVal] = calcMeanAndStd(raw)
    %   Compute mean and standard deviation for a cuboi of
    %   raw (i.e. grayscale) data.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    % convert to float
    raw = single(raw);

    % first, compute the mean value
    meanVal = mean(mean(mean(raw)));

    % now, the standard deviation
    % yes, this is a biased estimator
    sqDiff = (raw - meanVal) .^ 2;
    stdVal = sqrt(mean(mean(mean(sqDiff))));
end
