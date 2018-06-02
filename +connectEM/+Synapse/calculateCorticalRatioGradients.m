function [inhExcGradient, tcCcGradient] = ...
        calculateCorticalRatioGradients(param, conn, syn, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.binSize = 1;
    opt = Util.modifyStruct(opt, varargin{:});
    
    synT = connectEM.Connectome.buildSynapseTable(conn, syn);
    
    synPos = connectEM.Synapse.calculatePositions(param, syn);
    synPos = synPos - mean(param.bbox, 2)';
    synPos = synPos .* param.raw.voxelSize / 1E3;
    
    binEdges = ceil(max(abs(synPos(:))));
    binEdges = reshape((-binEdges):binEdges, [], 1);
    
    synData = [ ...
        discretize(synPos(synT.id, 1), binEdges), ...
        double(conn.axonMeta.axonClass(synT.preAggloId))];
    synData = accumarray(synData, 1, [numel(binEdges) - 1, 4]);
    
    inhExcGradient = calculateGradient(binEdges, synData, 3, 1:3);
    tcCcGradient = calculateGradient(binEdges, synData, 2, 1:2);
end

function fit = calculateGradient(binEdges, synData, numIds, denomIds)
    data = table;
    data.pos = ( ...
        binEdges(1:(end - 1)) ...
        + binEdges(2:end)) / 2;
    
    data.num = sum(synData(:, numIds), 2);
    data.denom = sum(synData(:, denomIds), 2);
    data.ratio = data.num ./ data.denom;
    
    weightIds = union(numIds, denomIds);
    data.weight = sum(synData(:, weightIds), 2);
    data(~data.weight, :) = [];
    
    A = [data.pos, ones(height(data), 1)];
    
    fitParam = (data.weight .* A) \ (data.weight .* data.ratio);
    fit = @(x) fitParam(1) .* x + fitParam(2);
end
