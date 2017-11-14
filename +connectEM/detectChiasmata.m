function output = detectChiasmata(p, nodesV, edges, outputFolder)
    % Detect chiasmata in skeletons based on marching sphere algorithm
    % Nodes should be in voxel, scaled here

    if ~isfield(p, 'minNrChiasmaExits') ...
            || isempty(p.minNrChiasmaExits)
        % for backward compatibility
        p.minNrChiasmaExits = 4;
    end

    % Create output folder if it does not exist
    if ~isempty(outputFolder) && ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Scale to nm
    % Make sure edges are unique
    nodes = bsxfun(@times, nodesV, p.voxelSize);
    edges = unique(edges, 'rows');

    % for each node ("marching sphere" approach to merger detection)
    nrExits = zeros(size(nodes, 1), 1);
    showProgress = @(n) fprintf( ...
        '%d of %d nodes done\n', n, size(nodes, 1));

    for i = 1:size(nodes, 1)
        nrExits(i) = connectEM.detectChiasmataNodes(p, nodes, edges, i);
        if ~mod(i, 10000); showProgress(i); end
    end

    % Mark intersections
    isIntersection = (nrExits >= p.minNrChiasmaExits);

    % Find CC of detected intersections according to graph
    output = ...
        connectEM.detectChiasmataPostSingleNodeLabel( ...
            p, nodes, nodesV, edges, isIntersection, nrExits);

    if ~isempty(outputFolder)
        fprintf('Writing result... ');
        Util.saveStruct(fullfile(outputFolder, 'result.mat'), output);
        fprintf('done!\n');
    end
end
