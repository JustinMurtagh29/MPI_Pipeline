function generateDendriteEndings(param,suffix)
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Christian Schramm <christian.schramm@brain.mpg.de>

    options = struct;
%     options.latentScore = 0.5;
    options.segDirScore = 0.8;
    options.distanceCutoff = 800; % in nm

    if ~exist('suffix','var')
        suffix = '';
    end
    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');

    % load directionality information
    endingInput = fullfile(dataDir, sprintf('dendriteEndingInputData%s.mat',suffix));
    endingInput = load(endingInput, 'dendriteIds', 'directionality');
    directionality = endingInput.directionality;
    dendriteIds = endingInput.dendriteIds;

    % load border CoMs
    borderCoM = fullfile(param.saveFolder, 'globalBorder.mat');
    borderCoM = load(borderCoM, 'borderCoM');
    borderCoM = borderCoM.borderCoM;

    % Find all borders for valid endings
%     idxDirectional = cellfun( ...
%         @(x) x(:, 1) > options.latentScore, ...
%         directionality.latent, 'UniformOutput', false);
    idxEnding = cellfun( ...
        @(x) abs(x) > options.segDirScore, ...
        directionality.scores, 'UniformOutput', false);
%     idxAll = cellfun( ...
%         @(x, y) find(x & y), ...
%         idxDirectional, idxEnding, 'UniformOutput', false);

    % Keep only those agglomerates with at least one ending candidate
    nrCanidates = cellfun(@numel, idxEnding);
    dendriteMask = nrCanidates > 0;
    dendriteIds = dendriteIds(dendriteMask);

    display([num2str(numel(dendriteMask)) ' agglomerates > 5 micron in total']);
    display([num2str(100 - (sum(dendriteMask)./numel(dendriteMask))*100, '%.2f') '% of > 5 micron agglomerates have not a single ending']);
    display([num2str(numel(dendriteMask) - sum(dendriteMask)) ' in total']);

    % clustering on left candidates
    borderIds = cellfun( ...
        @(x, y) x(y), ...
        directionality.borderIdx(dendriteMask), ...
        idxEnding(dendriteMask), 'UniformOutput', false);
    borderPositions = cellfun( ...
        @(x) borderCoM(x, :), borderIds, 'UniformOutput', false);
    borderClusters = cellfun( ...
        @(x) clusterBorders(param, options, x), ...
        borderPositions, 'UniformOutput', false);

    clusterSizes = cellfun(@max, borderClusters,'uni',0);
    singleEnding = sum(clusterSizes == 1);
    display([num2str(singleEnding./numel(clusterSizes)*100, '%.2f') '% of agglomerates have just one single ending']);
    display([num2str(singleEnding) ' in total']);

    % save all results for convenience in query generation
    gitInfo = Util.gitInfo();
    save(fullfile(dataDir, sprintf('dendriteEndingsAllData%s.mat',suffix)));

    % save result
    out = struct;
    out.dendriteIds = dendriteIds;
    out.dendriteMask = dendriteMask;
    out.borderIds = borderIds;
    out.borderPositions = borderPositions;
    out.borderClusters = borderClusters;
    out.gitInfo = Util.gitInfo();

    Util.saveStruct(fullfile(dataDir, sprintf('dendriteEndings%s.mat',suffix)), out);
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
