function job = buildSegmentMetaData(param)
    cubes = param.local;
    rootDir = param.saveFolder;
    
    jobSharedInputs = {param};
    jobInputs = arrayfun(@(curIdx) {{curIdx}}, 1:numel(cubes));

    % run job
    cluster = Cluster.getCluster('-l h_vmem=12G');
    job = Cluster.startJob( ...
        @buildInCube, jobInputs, ...
        'sharedInputs', jobSharedInputs, ...
        'cluster', cluster, 'name', mfilename());
    Cluster.waitForJob(job);
    
    % collect all local results
    loadMeta = @(p) load(fullfile(p.saveFolder, 'segmentMeta.mat'));
    meta = arrayfun(loadMeta, cubes, 'UniformOutput', false);
    meta = Util.concatStructs('last', meta{:});
    
    % find maximum segment ID
    meta.maxSegId = max(meta.maxSegId);
    
    % write global result
    metaFile = fullfile(rootDir, 'segmentMeta.mat');
    Util.saveStruct(metaFile, meta);
end

function buildInCube(param, cubeIdx)
    cubeParam = param.local(cubeIdx);
    cubeDir = cubeParam.saveFolder;
    
    % calculate meta data
    meta = Seg.Local.getSegmentMeta(param, cubeIdx);
    
    % save result
    metaFile = fullfile(cubeDir, 'segmentMeta.mat');
    Util.saveStruct(metaFile, meta);
end

