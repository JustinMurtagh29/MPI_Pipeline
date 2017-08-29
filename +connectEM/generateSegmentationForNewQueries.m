function generateSegmentationForNewQueries(param)

    % Define where to load and store results
    dataDir = fullfile(param.saveFolder, 'aggloState');

    % Load segmentMeta
    [~, segmentMeta] = connectEM.loadAllSegmentationData(param);

    % Load larger 5 micron agglomerates
    m = load(fullfile(dataDir, 'axons_04.mat'));
    axons = m.axons(m.indBigAxons);
    axons = arrayfun(@Agglo.fromSuperAgglo, axons, 'UniformOutput', false);
    clear m

    % Write new segmentation based on axon queries
    mapping = connectEM.createLookup(segmentMeta, axons);
    Seg.Global.applyMappingToSegmentation(param, mapping, fullfile(dataDir, '20170829_segmentationForQueries', '1'));

    % Create resolution pyramid
    thisBBox = [1, 1, 1; (ceil(param.bbox(:, 2) ./ 1024) .* 1024)']';
    seg.root = fullfile(param.saveFolder, 'aggloState', '20170829_segmentationForQueries', '1');
    seg.prefix = param.seg.prefix;
    createResolutionPyramid(seg, thisBBox, [], true);

end

