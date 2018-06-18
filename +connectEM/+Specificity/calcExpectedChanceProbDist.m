function [pVals, expAxonCount] = ...
        calcExpectedChanceProbDist(axonSynCounts, synProb, varargin)
    % Calculates the expected distribution of righ-tail p-values
    % (upper CDF), i.e. P(X >= p), for a binomial model with success rate
    % syn prob and axons with varying total synapse numbers.
    %
    % calcExpectedChanceProbDist(axonSynCounts, [], 'pdf', pdf_func)
    %   calculates the p-value distribution using the pdf given by
    %   pdf_func, where pdf_func is a function handle that accepts the
    %   total number of synapses n and must return a distribution over the
    %   0:n as a vector of probabilities.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    %   Benedikt Staffler <benedikt.staffler@brain.mpg.de>
    
    opts = struct;
    opts.pdf = @(n)binopdf(0:n, n, synProb);
    opts = Util.modifyStruct(opts, varargin{:});
    
   [synCounts, ~, axonCounts] = unique(axonSynCounts);
    axonCounts = accumarray(axonCounts, 1);
    
    pVals = cell(size(synCounts(:)));
    expAxonCount = cell(size(synCounts(:)));
    
    for curIdx = 1:numel(synCounts)
        curSynCount = synCounts(curIdx);
        curAxonCount = axonCounts(curIdx);
        
        curProbs = opts.pdf(curSynCount);
        curPVals = flip(cumsum(flip(curProbs)));
        
        pVals{curIdx} = curPVals(:);
        expAxonCount{curIdx} = curAxonCount * curProbs(:);
    end
    
    pVals = cell2mat(pVals);
    expAxonCount = cell2mat(expAxonCount);
    
   [pVals, ~, sortIds] = unique(pVals);
    expAxonCount = accumarray(sortIds, expAxonCount);
end
