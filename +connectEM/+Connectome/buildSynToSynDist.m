function [synToSynDists, synIds] = ...
        buildSynToSynDist(conn, syn, agglos, preOrPost, varargin)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    opt = struct;
    opt.segWeights = [];
    opt.showProgress = false;
    opt.voxelSize = ones(1, 3);
    opt = Util.modifyStruct(opt, varargin{:});
    
    switch preOrPost
        case 'pre'
            assert(isequal(size(conn.axons), size(agglos)));
            connAggloIds = conn.connectome.edges(:, 1);
            synAgglos = syn.synapses.presynId;
        case 'post'
            assert(isequal(size(conn.dendrites), size(agglos)));
            connAggloIds = conn.connectome.edges(:, 2);
            synAgglos = syn.synapses.postsynId;
    end
    
    if isempty(opt.segWeights)
        % By default, pick random segment
        opt.segWeights = max(cellfun(@max, synAgglos));
        opt.segWeights = ones(opt.segWeights, 1);
    end
    
    % Utility
    last = @(vals) vals(end);
    
    % Outputs
    synToSynDists = cell(size(agglos));
    synIds = cell(size(agglos));
    
    for curId = 1:numel(agglos)
        agglo = agglos(curId);
        
        % Find synapses
        curSynIds = connAggloIds == curId;
        curSynIds = conn.connectome.synIdx(curSynIds);
        curSynIds = cell2mat(curSynIds);
        
        % Map synapses onto node of largest segment
        curSynNodeIds = synAgglos(curSynIds);
        curSynNodeIds = cellfun( ...
            @(ids) Util.sortBy(ids, opt.segWeights(ids)), ...
            curSynNodeIds, 'UniformOutput', false);
        curSynNodeIds = cellfun(@(segIds) last(intersect( ...
            segIds, agglo.nodes(:, 4), 'stable')), curSynNodeIds);
       [~, curSynNodeIds] = ismember(curSynNodeIds, agglo.nodes(:, 4));
        
        % Calculate distances
        curSynToSynDist = SuperAgglo.pairwiseDist( ...
            agglo, curSynNodeIds, 'voxelSize', opt.voxelSize);
        
        synToSynDists{curId} = curSynToSynDist;
        synIds{curId} = curSynIds(:);
        
        if opt.showProgress; Util.progressBar(curId, numel(agglos)); end
    end
end
