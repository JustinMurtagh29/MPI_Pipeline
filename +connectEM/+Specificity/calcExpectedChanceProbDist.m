function [pVals, expAxonCount] = ...
        calcExpectedChanceProbDist(axonSynCounts, synProb)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
   [synCounts, ~, axonCounts] = unique(axonSynCounts);
    axonCounts = accumarray(axonCounts, 1);
    
    pVals = cell(size(synCounts(:)));
    expAxonCount = cell(size(synCounts(:)));
    
    for curIdx = 1:numel(synCounts)
        curSynCount = synCounts(curIdx);
        curAxonCount = axonCounts(curIdx);
        
        curProbs = binopdf(0:curSynCount, curSynCount, synProb);
        curPVals = flip(cumsum(flip(curProbs)));
        
        pVals{curIdx} = curPVals(:);
        expAxonCount{curIdx} = curAxonCount * curProbs(:);
    end
    
    pVals = cell2mat(pVals);
    expAxonCount = cell2mat(expAxonCount);
    
   [pVals, sortIds] = sort(pVals, 'ascend');
    expAxonCount = expAxonCount(sortIds);
end
