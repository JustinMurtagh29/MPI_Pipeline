function [probLower, probHigher] = ...
        calcRatioChanceProbs(classConn, numClassIds, denomClassIds)
    % calcRatioChanceProbs(classConn, numClassIds, denomClassIds)
    %   Given a class connectome and two sets of classes (A and B, with A
    %   being a subset of B), this function calculates for each neurite the
    %   probability of seeing a value A / B equally or more extreme than
    %   the observed one under assumption of a multinomial distribution.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Specificity.calcFractionProbs;
    
    % Sanity check
    assert(all(ismember(numClassIds, denomClassIds)));
    
    % Calculate observed fractions
    obsFracs = ...
        sum(classConn(:, numClassIds), 2) ...
     ./ sum(classConn(:, denomClassIds), 2);
    obsFracs(isnan(obsFracs)) = 0;
    
    % Calculate probabilities for multinomial
    otherClassIds = setdiff(1:size(classConn, 2), denomClassIds);
    denomClassIds = setdiff(denomClassIds, numClassIds);
    
    classProbs = sum(classConn, 1);
    classProbs = classProbs / sum(classProbs);
    
    mnProbs = [ ...
        sum(classProbs(numClassIds)), ...
        sum(classProbs(denomClassIds)), ...
        sum(classProbs(otherClassIds))];
    
    % Group by synapse count for efficiency
    synCounts = sum(classConn, 2);
   [uniSynCounts, ~, uniRows] = unique(synCounts);
   
    % Prepare output
    probLower = zeros(size(synCounts));
    probHigher = ones(size(synCounts));
    
    for curIdx = 1:numel(uniSynCounts)
        curRows = find(uniRows == curIdx);
        curSynCount = uniSynCounts(curIdx);
        
       [curFracs, curProbs] = ...
            calcFractionProbs(mnProbs, curSynCount);
        
        % Calculate probabilities
       [~, curThreshIdx] = ...
            ismember(obsFracs(curRows), curFracs);
        
        curProbs = cumsum(curProbs);
        probLower(curRows) = curProbs(curThreshIdx);
        
        curThreshIdx = curThreshIdx - 1;
        probHigher(curRows(curThreshIdx ~= 0)) = ...
            1 - curProbs(curThreshIdx(curThreshIdx ~= 0));
    end
end
