function [logAreas, gaussianIds] = sampleGmm(mix, varargin)
    % [logAreas, gaussianIds] = sampleGmm(mix, varargin)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opts = struct;
    opts.numSynapses = 2;
    opts.numConnections = 1;
    opts = Util.modifyStruct(opts, varargin{:});
    
    gaussianEdges = [mix.coeff];
    gaussianEdges = gaussianEdges / sum(gaussianEdges);
    gaussianEdges = cumsum([0; gaussianEdges(:)]);

    gaussianIds = linspace(0, 1, opts.numConnections);
    gaussianIds = discretize(gaussianIds, gaussianEdges);
    gaussianIds = gaussianIds(randperm(numel(gaussianIds)));
    gaussianIds = reshape(gaussianIds, [], 1);

    logAreas = randn(opts.numConnections, opts.numSynapses);
    logAreas = cat(1, mix(gaussianIds).std) .* logAreas;
    logAreas = cat(1, mix(gaussianIds).mean) + logAreas;
end
