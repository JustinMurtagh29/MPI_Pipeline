function pos = calculatePositions(param, syn, where)
    % pos = calculatePosition(param, syn, where)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('where', 'var') || isempty(where)
        where = 'prePost';
    end
    
    doRecenter = endsWith(where, 'Recenter');
    wmean = @(w, d) sum((w ./ sum(w)) .* d, 1);
    
    points = Seg.Global.getSegToCentroidMap(param);
    weights = Seg.Global.getSegToSizeMap(param);
    
    switch where
        case {'pre', 'preRecenter'}
            segIds = syn.synapses.presynId;
        case {'post', 'postRecenter'}
            segIds = syn.synapses.postsynId;
        case {'prePost', 'prePostRecenter'}
            segIds = cellfun( ...
                @vertcat, ...
                syn.synapses.presynId, ...
                syn.synapses.postsynId, ...
                'UniformOutput', false);
        otherwise
            error('Invalid input argument');
    end
    
    pos = cellfun( ...
        @(ids) wmean(weights(ids), points(ids, :)), ...
        segIds, 'UniformOutput', false);
    
    if doRecenter
        pos = cellfun( ...
            @(candPos, segIds) recenter( ...
                param, candPos, points(segIds, :)), ...
            pos, segIds, 'UniformOutput', false);
    end
    
    pos = cell2mat(pos);
end

function pos = recenter(param, candPos, segPos)
   [~, idx] = pdist2( ...
        segPos .* param.raw.voxelSize, ...
        candPos .* param.raw.voxelSize, ...
        'squaredeuclidean', 'Smallest', 1);
    
    assert(isscalar(idx));
    pos = segPos(idx, :);
end
