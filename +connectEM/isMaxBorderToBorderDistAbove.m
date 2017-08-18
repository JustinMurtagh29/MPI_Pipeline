function isAbove = isMaxBorderToBorderDistAbove(param, minDistNm, agglos)
    % isAbove = isMaxBorderToBorderDistAbove(param, minDistNm, agglos)
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    %% load data
    rootDir = param.saveFolder;
    graph = Graph.load(rootDir);
    
    % remove correspondences
    graph(~graph.borderIdx, :) = [];
    
    % find border for each agglo
    maxSegId = Seg.Global.getMaxSegId(param);
    aggloLUT = Agglo.buildLUT(maxSegId, agglos);
    
    graph.aggloIds = aggloLUT(graph.edges);
    graph(~diff(graph.aggloIds, 1, 2), :) = [];
    
    % group borders by agglomerate
    aggloBorderIds = accumarray( ...
        1 + graph.aggloIds(:), repmat(graph.borderIdx, 2, 1), ...
        [], @(borderIds) {borderIds});
    aggloBorderIds(1) = [];
    
    % load border CoM
    borderFile = fullfile(rootDir, 'globalBorder.mat');
    load(borderFile, 'borderCoM');
    
    % to physical units
    borderCoM = double(borderCoM); %#ok
    borderCoM = bsxfun(@times, borderCoM, param.raw.voxelSize);
    
    isAboveFunc = @(b) forBorders(borderCoM, minDistNm, b);
    isAbove = cellfun(isAboveFunc, aggloBorderIds);
end

function isAbove = forBorders(borderCoM, minDistNm, borderIds)
    isAbove = false;
    if numel(borderIds) < 2; return; end
    
    points = borderCoM(borderIds, :);
    
    pointsBox = [ ...
        min(points, [], 1); ...
        max(points, [], 1)];
    pointsDiff = diff(pointsBox);
    
    isAbove = any(pointsDiff > minDistNm);
    if isAbove; return; end;
    
    % do the precise calculation
    maxDist = pdist(points, 'squaredeuclidean');
    maxDist = sqrt(max(maxDist));
    
    isAbove = maxDist > minDistNm;
end
