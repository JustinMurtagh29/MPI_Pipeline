function output = detectChiasmata(p, nodesV, edges, outputFolder)
    % Detect chiasmata in skeletons based on marching sphere algorithm
    % Nodes should be in voxel, scaled here

    if ~isfield(p, 'minNrChiasmaExits') ...
            || isempty(p.minNrChiasmaExits)
        % for backward compatibility
        p.minNrChiasmaExits = 4;
    end
    
    if isempty(edges)
        % NOTE(amotta): It's entirely possible for an agglomerate to
        % consist of a single node and no edges. Unfortunately, somebody
        % messed up and forgot to make sure that `edges` has two columns.
        edges = reshape(edges, [], 2);
    end
    
    assert(size(nodesV, 2) == 3);
    assert(size(edges, 2) == 2);

    % Scale to nm
    % Make sure edges are unique
    nodes = bsxfun(@times, nodesV, p.raw.voxelSize);
    edges = unique(edges, 'rows');
    
    if size(nodes, 1) > 1E5
        % run locally
        nrExits = forLargeAgglo(p, nodes, edges);
    else
        % run in parallel on cluster
        nrExits = forSmallAgglo(p, nodes, edges);
    end

    % Mark intersections
    isIntersection = (nrExits >= p.minNrChiasmaExits);
    
    % Find CC of detected intersections according to graph
    output = ...
        connectEM.detectChiasmataPostSingleNodeLabel( ...
            p, nodes, nodesV, edges, isIntersection, nrExits);
    
    if ~isempty(outputFolder)
        if ~exist(outputFolder, 'dir')
            mkdir(outputFolder);
        end
        
        fprintf('Writing result... ');
        Util.save(fullfile(outputFolder, 'result.mat'), output);
        fprintf('done!\n');
    end
end

function nrExits = forSmallAgglo(p, nodes, edges)
    nrExits = forNodeIds(p, nodes, edges, 1:size(nodes, 1));
end

function nrExits = forLargeAgglo(p, nodes, edges)
    taskCount = 500;
    nodeCount = size(nodes, 1);
    blockSize = ceil(nodeCount / taskCount);
    
    inputArgs = arrayfun(@(off) {{ ...
        off:min(off + blockSize - 1, nodeCount)}}, ...
        1:blockSize:nodeCount);
    sharedArgs = {p, nodes, edges};
    
    %{
    cluster = Cluster.getCluster( ...
        '-p 0', ...
        '-pe openmp 1', ...
        '-l h_vmem=12G', ...
        '-l h_rt=2:00:00');
    %}

    cluster = Cluster.config( ...
        'scheduler', 'slurm', ...
        'memory', 12, ...
        'time', '2:00:00');
    
    job = Cluster.startJob( ...
        @forNodeIds, inputArgs, ...
        'sharedInputs', sharedArgs, ...
        'cluster', cluster, ...
        'numOutputs', 1);
    wait(job);
    
    nrExits = fetchOutputs(job);
    nrExits = sum(cat(2, nrExits{:, 1}), 2);
end

function nrExits = forNodeIds(p, nodes, edges, ids)
    nrExits = zeros(size(nodes, 1), 1);
    nrExits(ids) = arrayfun(@(i) ...
        connectEM.detectChiasmataNodes(p, nodes, edges, i), ids);
end
