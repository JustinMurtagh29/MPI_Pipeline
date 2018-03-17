function [expFrac, expCount] = ...
        calcExpectedFractionDist(classConn, numClassIds, denumClassIds)
    % calcExpectedFractionDist(classConn, numClassIds, denumClassIds)
    %   Given a class connectome and two sets of classes (A and B, with A
    %   being a subset of B), this function calculates the expected
    %   distribution of A / B.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Sanity check
    assert(all(ismember(numClassIds, denumClassIds)));
    
    otherClassIds = setdiff(1:size(classConn, 2), denumClassIds);
    denumClassIds = setdiff(denumClassIds, numClassIds);
    
    % Calculate probabilities for multinomial
    classProbs = sum(classConn, 1);
    classProbs = classProbs / sum(classProbs);
    
    mnProbs = [ ...
        sum(classProbs(numClassIds)), ...
        sum(classProbs(denumClassIds)), ...
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
        
        % Evaluate multinomial probability distribution.
        %
        % NOTE(amotta): The third entry becomes negative if the sum of
        % `curExc` and `curInh` is larger than the synapse count. This may
        % seem bad. But in fact it guarantees that `mnpdf` will return
        % zero and is thus very useful!
        curGrid = zeros(curSyns + 1, curSyns + 1, 3);
       [curGrid(:, :, 1), curGrid(:, :, 2)] = ndgrid(0:curSyns, 0:curSyns);
        curGrid(:, :, 3) = curSyns - sum(curGrid, 3);
        
        curGrid = reshape(curGrid, [], 3);
        curExpCount = mnpdf(curGrid, mnProbs);
        
        % Let's set e / (e + i) = 0 for e = i = 0.
        curExpFrac = curGrid(:, 1) ./ sum(curGrid(:, 1:2), 2);
        curExpFrac(isnan(curExpFrac)) = 0;
        
       [curExpFrac, ~, curUniIds] = unique(curExpFrac);
        curExpCount = accumarray(curUniIds, curExpCount);
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