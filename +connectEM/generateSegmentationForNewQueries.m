function generateSegmentationForNewQueries(param, seg, axFile)

    % Define where to load and store results
    dataDir = fullfile(param.saveFolder, 'aggloState');
    if ~exist('axFile', 'var') || isempty(axFile)
        axFile = fullfile(dataDir, 'axons_04.mat');
    end

    % Load segmentMeta
    [~, segmentMeta] = connectEM.loadAllSegmentationData(param);

    % Load larger 5 micron agglomerates
    if ischar(axFile)
        m = load(axFile);
        axons = m.axons(m.indBigAxons);
        axons = arrayfun(@Agglo.fromSuperAgglo, axons, 'uni', 0);
    else
        axons = axFile;
    end
    clear m;
    
    if ~exist('seg', 'var') || isempty(seg)
        seg = struct;
        seg.root = fullfile(param.saveFolder, 'aggloState', ...
            '20170829_segmentationForQueries', '1');
        seg.prefix = param.seg.prefix;
    end

    % Write new segmentation based on axon queries
    mapping = connectEM.createLookup(segmentMeta, axons);
    Seg.Global.applyMappingToSegmentation(param, mapping, seg);

    % Create resolution pyramid
    thisBBox = [1, 1, 1; (ceil(param.bbox(:, 2) ./ 1024) .* 1024)']';
    createResolutionPyramid(seg, thisBBox, [], true);
end
