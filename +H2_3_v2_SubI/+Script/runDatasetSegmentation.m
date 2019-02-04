function runDatasetSegmentation(p, pNew)
    
    info = Util.runInfo();
    
    %% Loading data
    meta = load(fullfile(p.saveFolder, 'segmentMeta.mat'));
    maxSegId = meta.maxSegId;
    points = transpose(meta.point);
    clear meta;
    
    % generate new segmentation per cube
    inputCell = arrayfun( ...
            @(newLocal, oldLocal) { ...
                newLocal.tempSegFile, [oldLocal.saveFolder 'graphH.mat']}, ...
            pNew.local(:), p.local(:), 'UniformOutput', false);
        
    job = Cluster.startJob( ...
            @jobWrapper, inputCell, ...
            'sharedInputs', {p, maxSegId}, ...
            'sharedInputsLocation', [1,2], ...
            'name', 'segNew', ...
            'cluster', { ...
            'memory', 36, ...
            'time', '24:00:00'});
end

function jobWrapper(p, maxSegId, saveFile, graphFile)
    Util.log('loading graph...')
    graph = load(graphFile));
    graph = graph.graph;

    % Agglomerate down to this score treshold
    minScore = 0;

    Util.log('Build agglomerates and export skeleton')
    mergeEdges = graph.edge(graph.score > minScore, :);
    [~, agglos] = Graph.buildConnectedComponents(maxSegId, mergeEdges);

    lut = Agglo.buildLUT(maxSegId, agglos);
    lut(lut == 0) = max(lut) + (1:sum(lut == 0));
    lut = [0; lut(:)];
    segNew = lut(seg + 1);
    Util.save(saveFile, segNew);

end
