function boutonAgglos = buildBoutonAgglos( ...
        segPoints, conn, syn, interSyn, varargin)
    import connectEM.Axon.*;
    
    opt = struct;
    opt.voxelSize = [];
    opt.distThresh = 1000;
    opt.parallelize = false;
    opt.showProgressBar = false;
    opt = Util.modifyStruct(opt, varargin{:});
    
    if ~isempty(opt.voxelSize)
        assert(isequal(size(opt.voxelSize), [1, 3]));
        segPoints = opt.voxelSize .* segPoints;
    end
    
    % Build boutons via synapses
    synIds = getSynapses(conn, syn);
    boutonIds = clusterSynapsesIntoBoutons(synIds, interSyn);
    
    % Find segments
    axonAgglos = conn.axons;
    preSynAgglos = syn.synapses.presynId;
    boutonAgglos = cell(numel(axonAgglos), 1);
    
    if opt.parallelize
        parfor curId = 1:numel(axonAgglos)
            boutonAgglos{curId} = forAxon( ...
                opt.distThresh, segPoints, axonAgglos{curId}, ...
                preSynAgglos(synIds{curId}), boutonIds{curId}); %#ok
        end
    else
        tic;
        for curId = 1:numel(axonAgglos)
            boutonAgglos{curId} = forAxon( ...
                opt.distThresh, segPoints, axonAgglos{curId}, ...
                preSynAgglos(synIds{curId}), boutonIds{curId});
            
            if ~opt.showProgressBar; continue; end
            Util.progressBar(curId, numel(axonAgglos));
        end
    end
end

function boutonAgglos = forAxon( ...
        distThresh, segPoints, segIds, preSynAgglos, boutonIds)
    % Graph representation for distance along neurite
    skel = SuperAgglo.fromAgglo( ...
        {segIds}, segPoints, 'mst');
    
    edgeLengths = ...
        skel.nodes(skel.edges(:, 1), 1:3) ...
      - skel.nodes(skel.edges(:, 2), 1:3);
    edgeLengths = sqrt(sum(edgeLengths .^ 2, 2));
    
    skel = graph( ...
        skel.edges(:, 1), skel.edges(:, 2), ...
        edgeLengths, size(skel.nodes, 1));
    
    % KD-tree for first-order limit on bouton size
    kdTree = KDTreeSearcher(segPoints(segIds, :));

    boutonCount = max(boutonIds);
    boutonAgglos = cell(boutonCount, 1);

    for curBoutonId = 1:boutonCount
        curSynSegIds = preSynAgglos(boutonIds == curBoutonId);
        curSynSegIds = unique(cell2mat(curSynSegIds));

        curSegIds = rangesearch( ...
            kdTree, segPoints(curSynSegIds, :), distThresh);
        curSegIds = cell2mat(reshape(curSegIds, 1, []));
        curSegIds = unique(curSegIds(:));

       [~, curSynSegIds] = ismember(curSynSegIds, segIds);
        curSynSegIds = setdiff(curSynSegIds, 0);

        curDists = distances(skel, curSynSegIds, curSegIds);
        curSegIds = curSegIds(any(curDists < distThresh, 1));
        boutonAgglos{curBoutonId} = segIds(curSegIds(:));
    end
end
