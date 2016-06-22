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
    
    % skip borders across cubes
    dropMask = ~cubeIds;
    cubeIds(dropMask) = [];
    edges(dropMask, :) = [];
    
    % make continuous-relabeling
    [uniCubeIds, ~, cubeIds] = unique(cubeIds);
    
    % task params
    uniCubeCount = numel(uniCubeIds);
    taskParams = cell(uniCubeCount, 1);
    
    % load features
    for curIdx = 1:uniCubeCount
        curCubeIdx = uniCubeIds(curIdx);
        curCubeParam = param.local(curCubeIdx);
        
        % get edges
        curEdgeMask = (cubeIds == curIdx);
        curEdges = edges(curEdgeMask, :);
        
        % build parameters
        taskParams{curIdx} = { ...
            curCubeParam, curEdges};
    end
    
    % create job
    global CLUSTER_CPU;
    job = createJob(CLUSTER_CPU);
    job.Name = 'loadFeatures';
    
    % start and wait for job
    createTask(job, @loadFeaturesCube, 3, taskParams);
    
    % wait for job
    submit(job);
    wait(job);
    
    % get output
    out = fetchOutputs(job);
    
    % build output
    outEdges = vertcat(out{:, 1});
    outFeats = vertcat(out{:, 2});
    outFeatNames = out{end, 3};
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