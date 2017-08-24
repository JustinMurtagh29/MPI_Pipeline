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
    
    % Keep only agglomerates with at least one ending candidate
    nrCanidates = cellfun(@numel, idxAll);
    axonMask = nrCanidates > 0;
    axonIds = find(axonMask);
    display([num2str(numel(axonMask)) ' agglomerates in total']);
    display([num2str(1 - (sum(axonMask)./numel(axonMask))*100, '%.2f') '% of agglomerates have not a single ending']);
    display([num2str(numel(axonMask) - sum(axonMask)) ' in total']);
    
        
    % clustering on left candidates
    borderIds = cellfun( ...
        @(x, y) x(y), ...
        directionality.borderIdx(axonMask), ...
        idxAll(axonMask), 'UniformOutput', false);
    borderPositions = cellfun( ...
        @(x) borderCoM(x, :), borderIds, 'UniformOutput', false);
    borderClusters = cellfun( ...
        @(x) clusterBorders(param, options, x), ...
        borderPositions, 'UniformOutput', false);
    
    clusterSizes = cellfun(@max, borderClusters);
    singleEnding = sum(clusterSizes == 1);
    display([num2str(singleEnding./numel(clusterSizes)*100, '%.2f') '% of agglomerates have just one single ending']);
    display([num2str(singleEnding) ' in total']);
    
    
    % save result
    out = struct;
    out.axonMask = axonMask;
    out.axonIds = axonIds;
    out.idxAll = idxAll;
    out.borderIds = borderIds;
    out.borderPositions = borderPositions;
    out.borderClusters = borderClusters;
    out.gitInfo = Util.gitInfo();
    
    Util.saveStruct(fullfile(param.saveFolder, 'aggloState/', 'axonEndings.mat'), out);
end

function clusterIds = clusterBorders(param, options, borderCoM)
    if size(borderCoM, 1) > 1
        % convert borderCoM to nm space
        borderCoM = double(borderCoM);
        borderCoM = bsxfun(@times, borderCoM, param.raw.voxelSize);
        
        clusterIds = clusterdata( ...
            borderCoM, ...
            'linkage', 'single', ...
            'criterion', 'distance', ...
            'cutoff', options.distanceCutoff);
    else
        clusterIds = ones(size(borderCoM, 1), 1);
    end
end
