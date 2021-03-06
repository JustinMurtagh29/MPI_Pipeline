function [expFrac, expCount] = ...
        calcExpectedRatioDist(classConn, numClassIds, denomClassIds)
    % calcExpectedRatioDist(classConn, numClassIds, denumClassIds)
    %   Given a class connectome and two sets of classes (A and B, with A
    %   being a subset of B), this function calculates the expected
    %   distribution of A / B.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    import connectEM.Specificity.calcRatioProbs;
    
    % Sanity check
    assert(all(ismember(numClassIds, denomClassIds)));
    
    otherClassIds = setdiff(1:size(classConn, 2), denomClassIds);
    denomClassIds = setdiff(denomClassIds, numClassIds);
    
    % Calculate probabilities for multinomial
    classProbs = sum(classConn, 1);
    classProbs = classProbs / sum(classProbs);
    
    mnProbs = [ ...
        sum(classProbs(numClassIds)), ...
        sum(classProbs(denomClassIds)), ...
        sum(classProbs(otherClassIds))];
    
    % Determine minimal set of synapse counts
    synCounts = sum(classConn, 2);
   [uniSyns, ~, uniSynNeurites] = unique(synCounts);
    uniSynNeurites = accumarray(uniSynNeurites, 1);
    
    expFrac = cell(numel(uniSyns), 1);
    expCount = cell(numel(uniSyns), 1);
    
    for curIdx = 1:numel(uniSyns)
        curSyns = uniSyns(curIdx);
        curNeurites = uniSynNeurites(curIdx);
        
       [curExpFrac, curExpCount] = calcRatioProbs(mnProbs, curSyns);
        curExpCount = curNeurites * curExpCount;
        
        expFrac{curIdx} = curExpFrac;
        expCount{curIdx} = curExpCount;
    end
    
    % Build output
    expFrac = cell2mat(expFrac);
   [expFrac, ~, uniIds] = unique(expFrac);
   
    expCount = cell2mat(expCount);
    expCount = accumarray(uniIds, expCount);
end
