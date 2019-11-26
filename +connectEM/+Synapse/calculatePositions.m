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
    
    coms = Seg.Global.getSegToCentroidMap(param);
    weights = Seg.Global.getSegToSizeMap(param);
    points = Seg.Global.getSegToPointMap(param);
    
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
        case 'border'
            edgeToBorder = fullfile(param.saveFolder, 'graph.mat');
            edgeToBorder = Util.load(edgeToBorder, 'borderIdx');
            
            segIds = cellfun( ...
                @(edgeIds) edgeToBorder(edgeIds), ...
                syn.synapses.edgeIdx, 'UniformOutput', false);
            
            coms = fullfile(param.saveFolder, 'globalBorder.mat');
           [coms, weights] = Util.load(coms, 'borderCoM', 'borderSize');
           
            coms = double(coms);
            weights = double(weights);
            points = coms;
        otherwise
            error('Invalid input argument');
    end
    
    pos = nan(numel(segIds), 3);
    for curIdx = 1:numel(segIds)
        curSegIds = segIds{curIdx};
        curPos = wmean(weights(curSegIds), coms(curSegIds, :));
        
        if doRecenter
            % Move from center of mass to closest segment position
            curPos = recenter(param, curPos, points(curSegIds, :));
        end
        
        pos(curIdx, :) = curPos;
    end
end

function pos = recenter(param, pos, segPos)
    if isequal(pos, segPos); return; end
    pos = (segPos - pos) .* param.raw.voxelSize;
   [~, pos] = min(sum(pos .* pos, 2));
    pos = segPos(pos, :);
end
