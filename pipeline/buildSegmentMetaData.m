function buildSegmentMetaData(param)
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    cubes = param.local;
    rootDir = param.saveFolder;
    
    jobSharedInputs = {param};
    jobInputs = arrayfun( ...
        @(curIdx) {curIdx}, 1:numel(cubes), ...
        'UniformOutput', false);

    % run job
    job = Cluster.startJob( ...
        @buildInCube, jobInputs, ...
        'sharedInputs', jobSharedInputs, ...
        'name', mfilename(), ...
        'cluster', { ...
            'memory', 12, ...
            'time', '00:30:00'});
    Cluster.waitForJob(job);
    
    % collect all local results
    Util.log('Building global segment meta data file');
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

