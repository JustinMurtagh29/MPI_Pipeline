function features = collectSegFeatures(param, meta)
    % segIds = collectSegFeatures(param, coords)
    %   This function collects the feature data of segIds given in
    %   meta.segIds
    %
    % param
    %   Parameter structure
    %
    % meta
    %   meta structure for which the features are wanted
    %
    % features
    %   NxM vector. The entry features(i,:) contains the features of the
    %   segment ID of meta.segIds(i).
    %
    % Written by
    %   Marcel Beining <marcel.beining@brain.mpg.de>

    
    load(fullfile(param.saveFolder,'SegmentFeatureMap.mat'),'fm')
    
    % build cube ids
    [uniCubeIds, ~, Idx2uniCube] = unique(meta.cubeIdx);
    
    % look up cubes
    uniCubeCount = size(uniCubeIds, 1);
    
    % prepare output
    features = NaN(numel(meta.segIds), fm.numFeatures);
    
    tic;
    for curIdx = 1:uniCubeCount
        curIdx2uniCube = Idx2uniCube == curIdx;
        % load features
        thisfeat = load(fullfile(param.local(uniCubeIds(curIdx)).saveFolder,'SegmentFeatures.mat'));
        thisMeta = load(fullfile(param.local(uniCubeIds(curIdx)).saveFolder,'segmentMeta.mat'));
        % find the position of the segIds which have volume above areaT
        [~,ia] = intersect(thisMeta.segIds(thisMeta.voxelCount>fm.areaT),meta.segIds(curIdx2uniCube));
        
        % delete indices to meta which are below area threshold
        curIdx2uniCube(meta.voxelCount <= fm.areaT) = 0; 
        % put the features of the segIds into the feature vector
        features(curIdx2uniCube,:) = thisfeat.features(ia,:);
        
        Util.progressBar(curIdx, uniCubeCount);
    end
end