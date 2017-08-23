function generateAxonEndings(param)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    
    options = struct;
    options.latentScore = 0.5;
    options.segDirScore = 0.8;
    options.distanceCutoff = 600; % in nm
    
    % load directionality information
    directionality = fullfile(param.saveFolder, 'axonEndingInputData.mat');
    directionality = load(directionality, 'directionality');
    directionality = directionality.directionality;
    
    % load border CoMs
    borderCoM = fullfile(param.saveFolder, 'globalBorder.mat');
    borderCoM = load(borderCoM, 'borderCoM');
    borderCoM = borderCoM.borderCoM;
    
    % convert border CoMs to nm space
    borderCoM = double(borderCoM);
    borderCoM = bsxfun(@times, borderCoM, param.raw.voxelSize);
    
    % Find all borders for valid endings
    idxDirectional = cellfun( ...
        @(x) x(:, 1) > options.latentScore, ...
        directionality.latent, 'UniformOutput', false);
    idxEnding = cellfun( ...
        @(x) abs(x) > options.segDirScore, ...
        directionality.scores, 'UniformOutput', false);
    idxAll = cellfun( ...
        @(x, y) find(x & y), ...
        idxDirectional, idxEnding, 'UniformOutput', false);
    
    % Keep only largest score in each agglomerate for now
    nrCanidates = cellfun(@numel, idxAll);
    axonMask = nrCanidates > 0;
    axonIds = find(axonMask);
    
    % clustering on left candidates
    borderIds = cellfun( ...
        @(x, y) x(y), ...
        directionality.borderIdx(axonMask), ...
        idxAll(axonMask), 'UniformOutput', false);
    borderPositions = cellfun( ...
        @(x) borderCoM(x, :), borderIds, 'UniformOutput', false);
    borderClusters = cellfun( ...
        @(x) clusterBorders(options, x), ...
        borderPositions, 'UniformOutput', false);
    
    % save result
    out = struct;
    out.axonMask = axonMask;
    out.axonIds = axonIds;
    out.borderIds = borderIds;
    out.borderPositions = borderPositions;
    out.borderClusters = borderClusters;
    
    Util.saveStruct(fullfile(param.saveFolder, 'axonEndings.mat'), out);
end

function clusterIds = clusterBorders(options, borderCoM)
    if size(borderCoM, 1) > 1
        clusterIds = clusterdata( ...
            borderCoM, ...
            'linkage', 'single', ...
            'criterion', 'distance', ...
            'cutoff', options.distanceCutoff);
    else
        clusterIds = ones(size(borderCoM, 1), 1);
    end
end
