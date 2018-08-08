function pos = calculatePositions(param, syn, where)
    % pos = calculatePosition(param, syn, where)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    if ~exist('where', 'var') || isempty(where)
        where = 'prePost';
    end
    
    points = Seg.Global.getSegToCentroidMap(param);
    weights = Seg.Global.getSegToSizeMap(param);
    
    switch where
        case 'pre'
            pos = syn.synapses.presynId;
        case 'post'
            pos = syn.synapses.postsynId;
        case 'prePost'    
            pos = cellfun( ...
                @vertcat, ...
                syn.synapses.presynId, ...
                syn.synapses.postsynId, ...
                'UniformOutput', false);
        otherwise
            error('Invalid input argument');
    end

    wmean = @(w, d) sum((w ./ sum(w)) .* d, 2);
    pos = cell2mat(cellfun( ...
        @(ids) wmean(weights(ids), points(ids, :)), ...
        pos, 'UniformOutput', false));
end
