function [outEdges, outFeats, outFeatNames] = ...
        loadFeatures(param, edges)
    % loadFeatures(param, edges)
    %   Loads the features for all of the specified edges.
    %
    % param
    %   Parameter structure produced by setParameterSettings
    %
    % edges
    %   Nx2 matrix. Each rows contains a pair of global seg-
    %   ment IDs that make up an edge.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % load seg-to-cube map
    rootDir = param.saveFolder;
    segToCubeFile = [rootDir, 'segToCubeMap.mat'];
    load(segToCubeFile, 'segToCubeMap');
    
    % find cube ids
    cubeIds = findCubeIds(segToCubeMap, edges);
    [uniCubeIds, ~, cubeIds] = unique(cubeIds);
    
    % skip borders across cubes
    if ~uniCubeIds(1)
        uniCubeIds = uniCubeIds(2:end);
        cubeIds = cubeIds - 1;
    end
    
    % prepare output
    uniCubeCount = numel(uniCubeIds);
    outEdges = cell(uniCubeCount, 1);
    outFeats = cell(uniCubeCount, 1);
    outFeatNames = cell(0);
    
    % load features
    tic;
    for curIdx = 1:uniCubeCount
        curCubeIdx = uniCubeIds(curIdx);
        curCubeParam = param.local(curCubeIdx);
        
        % get edges
        curEdgeMask = (cubeIds == curIdx);
        curEdges = edges(curEdgeMask, :);
        
        % handle current cube
        [curEdges, curFeats, curFeatNames] = ...
            loadFeaturesCube(curCubeParam, curEdges);
        
        % write back result
        outEdges{curIdx} = curEdges;
        outFeats{curIdx} = curFeats;
        outFeatNames = curFeatNames;
        
        % show progress
        Util.progressBar(curIdx, uniCubeCount);
    end
    
    % build output
    outEdges = vertcat(outEdges{:});
    outFeats = vertcat(outFeats{:});
end

function [outEdges, outFeats, outFeatNames] = ...
        loadFeaturesCube(cubeParam, edges)
    cubeDir = cubeParam.saveFolder;
    
    edgeFile = [cubeDir, 'edges.mat'];
    borderFeatFile = [cubeDir, 'bordersExtFeats.mat'];
    
    % load edges and border features
    cubeEdges = load(edgeFile, 'edges');
    cubeFeats = load(borderFeatFile);
    
    % build mask
    keepMask = ismember( ...
        cubeEdges.edges, edges, 'rows');
    
    % build output
    outEdges = cubeEdges.edges(keepMask, :);
    outFeats = cubeFeats.feats(keepMask, :);
    outFeatNames = cubeFeats.featNames;
end

function out = findCubeIds(segToCubeMap, edges)
    % apply mapping to each segment
    cubeIds = segToCubeMap(edges);
    
    % keep edges within same cube
    sameMask = (cubeIds(:, 1) == cubeIds(:, 2));
    
    % prepare output
    edgeCount = size(edges, 1);
    
    % build output
    out = zeros(edgeCount, 1);
    out(sameMask) = cubeIds(sameMask, 1);
end