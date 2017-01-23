function job = buildSegmentMetaData(param)
    cubes = param.local;
    rootDir = param.saveFolder;

    taskCount = numel(cubes);
    jobParams = arrayfun(@(curIdx) {{param, curIdx}}, 1:taskCount);
    jobParams = cellfun(@(x){x}, jobParams, 'uni', 0);

    % run job
    job = startCPU(@buildInCube, jobParams, mfilename());
    Cluster.waitForJob(job);
    
    % collect all local results
    loadMeta = @(p) load(fullfile(p.saveFolder, 'segmentMeta.mat'));
    meta = arrayfun(loadMeta, cubes, 'UniformOutput', false);
    meta = Util.concatStructs('last', meta{:});
    
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
