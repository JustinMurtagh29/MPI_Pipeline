function firstHitProbs = calcFirstHitProbs(classConn)
    % firstHitProbs = calcFirstHitProbs(classConn)
    %   Calculate the probability `firstHitProbs(i)` which best explains
    %   wheter or not the target class in `classConn(:, i)` is innervated
    %   under a binomial distribution.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    syns = sum(classConn, 2);

    optimizationOptions = optimoptions( ...
        'lsqnonlin', 'display', 'off', ...
        'SpecifyObjectiveGradient', true);
    
    firstHitProbs = nan(1, size(classConn, 2));
    for curIdx = 1:numel(firstHitProbs)
        curHits = classConn(:, curTargetId) > 0;

        curFuncs = @(x) (1 - curHits) + (2 .* curHits - 1) .* (x .^ syns);
        curDerivs = @(x) (2 .* curHits - 1) .* syns .* (x .^ (syns - 1));
        
        curPInv = lsqnonlin( ...
            @(x) deal(curFuncs(x), curDerivs(x)), ...
            0.5, 0, 1, optimizationOptions);
        
        firstHitProbs(curIdx) = 1 - curPInv;
    end
end

