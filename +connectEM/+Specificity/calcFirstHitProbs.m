function firstHitProbs = calcFirstHitProbs(classConn, method)
    % firstHitProbs = calcFirstHitProbs(classConn)
    %   Calculate the probability `firstHitProbs(i)` which best explains
    %   wheter or not the target class in `classConn(:, i)` is innervated
    %   under a binomial distribution.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('method', 'var') || isempty(method)
        % NOTE(amotta): For backward compatibility
        method = 'oneVersusRestBinomial';
    end
    
    switch method
        case 'oneVersusRestBinomial'
            firstHitProbs = oneVersusRestBinomial(classConn);
        case 'multinomial'
            firstHitProbs = multinomial(classConn);
        otherwise
            error('Unknown method "%s"', method);
    end
    
    firstHitProbs = reshape(firstHitProbs, 1, []);
end

function firstHitProbs = oneVersusRestBinomial(classConn)
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

function firstHitProbs = multinomial(classConn)
    % NOTE(amotta): Initial guess of first-hit probability is that all
    % target classes are innervated with equal probability.
    syns = sum(classConn, 2);
    classConn = classConn > 0;
    
    N = size(classConn, 2);
    x0 = repelem(1 / N, N, 1);
    
    % NOTE(amotta): This inequality-constraint makes sure that all first-
    % hit probabilities are non-negative.
    A = -eye(N);
    B = zeros(N, 1);
    
    % NOTE(amotta): This equality-constrait ensures that the sum of all
    % first-hit probabilities adds up to one. In combination with the above
    % inequality-constraint this also guarantees that the individual first-
    % hit probabilities are constrained to the interval from zero to one.
    Aeq = ones(1, N);
    Beq = 1;
    
    F = @(p) -sum(sum( ...
        log(1 - (1 - p') .^ syns) .* classConn ...
      + log((1 - p') .^ syns) .* (1 - classConn)));
  
    opts = optimoptions(@fmincon, 'Display', 'notify');
    firstHitProbs = fmincon(F, x0, A, B, Aeq, Beq, [], [], [], opts);
end
