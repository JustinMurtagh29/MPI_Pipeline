function firstHitProbs = calcFirstHitProbs(classConn)
    % firstHitProbs = calcFirstHitProbs(classConn)
    %   Calculate the probability `firstHitProbs(i)` which best explains
    %   wheter or not the target class in `classConn(:, i)` is innervated
    %   under a binomial distribution.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    syns = sum(classConn, 2);
    
    firstHitProbs = nan(1, size(classConn, 2));
    for curTargetId = 1:numel(firstHitProbs)
        curHits = classConn(:, curTargetId) > 0;
        
        curFunc = @(pTilde) -sum( ...
            (1 - curHits) .* log(pTilde .^ syns) ...
            + curHits .* log(1 - pTilde .^ syns));

        curPTilde = fminbnd(curFunc, 0, 1);
        firstHitProbs(curTargetId) = 1 - curPTilde;
    end
end

