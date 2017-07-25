function job = buildSegmentMetaData(param)
    cubes = param.local;
    rootDir = param.saveFolder;
    
    jobSharedInputs = {param};
    jobInputs = arrayfun(@(curIdx) {{curIdx}}, 1:numel(cubes));

    % run job
    cluster = Cluster.getCluster('-l h_vmem=12G -l s_rt=00:29:30 -l h_rt=00:29:55');
    job = Cluster.startJob( ...
        @buildInCube, jobInputs, ...
        'sharedInputs', jobSharedInputs, ...
        'cluster', cluster, 'name', mfilename());
    Cluster.waitForJob(job);
    
    % collect all local results
    loadMeta = @(p) load(fullfile(p.saveFolder, 'segmentMeta.mat'));
    meta = arrayfun(loadMeta, cubes, 'UniformOutput', false);
    meta = Util.concatStructs('last', meta{:});
%     meta = Util.concatStructs(struct('segIds',1,'voxelCount',1,'box',3,'centroid',2,'point',2,'maxSegId',1,'atborder',1,'cubeIdx',1) ,meta{:});
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

