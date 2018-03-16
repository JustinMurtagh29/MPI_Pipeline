function [expFrac, expCount] = ...
        calcExpectedFractionDist(classConn, numClassId, otherClassId)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % Calculate probabilities for multinomial
    classProbs = sum(classConn, 1);
    classProbs = classProbs / sum(classProbs);
    
    mnProbs = classProbs([numClassId, otherClassId]);
    mnProbs = horzcat(mnProbs, 1 - sum(mnProbs));
    
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