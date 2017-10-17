function generateSegmentationForNewDendriteQueries(param)

    % Define where to load and store results
    dataDir = fullfile(param.saveFolder, 'aggloState');
    segDir = '/tmpscratch/scchr';

    % Load segmentMeta
    [~, segmentMeta] = connectEM.loadAllSegmentationData(param);

    % Load larger 5 micron agglomerates
    m = load(fullfile(dataDir, 'dendrites_03_v2_splitmerged.mat'));
    dendrites = m.dendrites(m.indBigDends);
    dendrites = arrayfun(@Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);
    clear m

    % Write new segmentation based on axon queries
    mapping = connectEM.createLookup(segmentMeta, dendrites);
    Seg.Global.applyMappingToSegmentation(param, mapping, fullfile(segDir, '20171017_segmentationForDendriteQueries', '1'));

    % Create resolution pyramid
    thisBBox = [1, 1, 1; (ceil(param.bbox(:, 2) ./ 1024) .* 1024)']';
    seg.root = fullfile(segDir, '20171017_segmentationForDendriteQueries', '1');
    seg.prefix = param.seg.prefix;
    createResolutionPyramid(seg, thisBBox, [], true);

end

