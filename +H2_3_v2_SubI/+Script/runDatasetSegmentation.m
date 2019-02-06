% load old p and new p parameters and save new segmentation after Hierarchichal clustering
Util.log('Loading old and new pipeline parameter...')
m = load('/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet/allParameter.mat');
p = m.p;

m = load('/tmpscratch/sahilloo/data/H2_3_v2_U1_SubI/pipelineRun_mr2e_wsmrnet_HC_v2/allParameter.mat');
pNew = m.p;

Util.log('submitting job for segmentation writing per cube')
job = runSegmentation(p, pNew);

% wait for job to finish
Cluster.waitForJob(job,'',true);

Util.log('run pipeline steps for segmentation generation')
runPipeline(pNew, PipelineStep.OverlapRemoval, PipelineStep.CompressSegmentation)

function job = runSegmentation(p, pNew)
    
    info = Util.runInfo();
    
    %% Loading data
    meta = load(fullfile(p.saveFolder, 'segmentMeta.mat'));
    maxSegId = meta.maxSegId;
    clear meta;

    margin = p.tileBorder(:, 2);
    assert(isequal(margin, [256; 256; 128]));

    % generate new segmentation per cube
    inputCell = arrayfun( ...
            @(newLocal, oldLocal) { ...
                newLocal.tempSegFile, [oldLocal.saveFolder 'graphH.mat'] ...
                oldLocal.bboxSmall}, ...
            pNew.local(:), p.local(:), 'UniformOutput', false);
        
    job = Cluster.startJob( ...
            @jobWrapper, inputCell, ...
            'sharedInputs', {p, maxSegId, margin}, ...
            'sharedInputsLocation', [1,2,3], ...
            'name', 'segNew', ...
            'cluster', { ...
            'memory', 36, ...
            'time', '24:00:00'});
end

function jobWrapper(p, maxSegId, margin, saveFile, graphFile, bboxSmall)
    Util.log('loading graph...')
    graph = load(graphFile);
    graph = graph.graph;

    Util.log('loading segmentation...')
    boxLarge = bboxSmall + [-1, +1] .* margin(:);
    boxLarge = max(boxLarge, p.bbox(:, 1));
    boxLarge = min(boxLarge, p.bbox(:, 2));
    segOld = loadSegDataGlobal(p.seg, boxLarge);

    % Agglomerate down to this score treshold
    minScore = 0;

    Util.log('Build agglomerates and export skeleton')
    mergeEdges = graph.edge(graph.score > minScore, :);
    [~, agglos] = Graph.buildConnectedComponents(maxSegId, mergeEdges);

    lut = Agglo.buildLUT(maxSegId, agglos);
    lut(lut == 0) = max(lut) + (1:sum(lut == 0));
    lut = [0; lut(:)];
    seg = lut(segOld + 1);

    saveFolder = fileparts(saveFile);
    if ~exist(saveFolder, 'dir')
        mkdir(saveFolder);
    end
    Util.save(saveFile, seg);

end
