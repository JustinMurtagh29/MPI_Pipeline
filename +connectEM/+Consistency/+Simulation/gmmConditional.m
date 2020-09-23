function [sizeCondSizeProbs, sizeCondMixProbs, sizeProbs] = ...
        gmmConditional(mix, varargin)
    % condMap = gmmConditional(mix, varargin)
    %   This function generates a map of the conditional size distribution
    %   for same-axon same-dendrite synapse pairs which are simulated by a
    %   mixture of Gaussian distributions.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.xVec = [];
    opts.yVec = [];
    opts = Util.modifyStruct(opts, varargin{:});
    
    if xor(isempty(opts.xVec), isempty(opts.yVec))
        opts.xVec = union(opts.xVec, opts.yVec);
        opts.yVec = union(opts.xVec, opts.yVec);
    end
    
    opts.xVec = reshape(sort(opts.xVec), 1, []);
    opts.yVec = reshape(sort(opts.yVec), 1, []);
    
    sq = @(x) x .* x;
    mu = cat(3, mix.mean);
    sigmaSq = sq(cat(3, mix.std));
    p = cat(3, mix.coeff);
    p = p / sum(p);
    
    sizeCondMixProbs = ...
        p(:) ./ sqrt(2 * pi * sigmaSq(:)) ...
     .* exp(-sq(opts.xVec - mu(:)) ./ (2 * sigmaSq(:)));
 
    sizeProbs = sum(sizeCondMixProbs, 1);
    sizeCondMixProbs = sizeCondMixProbs ./ sizeProbs;
    
    sizeCondSizeProbs = permute( ...
        sizeCondMixProbs, [3, 2, 1]);
    sizeCondSizeProbs = ...
        sizeCondSizeProbs ./ sqrt(2 * pi * sigmaSq) ...
     .* exp(-sq(opts.yVec(:) - mu) ./ (2 * sigmaSq));
    sizeCondSizeProbs = sum(sizeCondSizeProbs, 3);
end
