function generateSegmentationForNewDendriteQueries(param)

    % Define where to load and store results
    dataDir = fullfile(param.saveFolder, 'aggloState');
    segDir = '/tmpscratch/scchr/20171027_segmentationForDendriteQueries';

    % Load segmentMeta
    [~, segmentMeta] = connectEM.loadAllSegmentationData(param);

    % Load larger 5 micron agglomerates
    m = load(fullfile(dataDir, 'dendrites_flight_01.mat'));
    dendrites = m.dendrites(m.indBigDends);
    dendrites = arrayfun(@Agglo.fromSuperAgglo, dendrites, 'UniformOutput', false);
    mapping = connectEM.createLookup(segmentMeta, dendrites);
    clear m
    
    seg = struct;
    seg.root = fullfile(segDir, '1');
    seg.prefix = param.seg.prefix;
    seg.backend = 'wkwrap';
    
    % Initialize WKW dataset
    wkwInit('new', seg.root, 32, 32, 'uint32', 1);
    Seg.Global.applyMappingToSegmentation(param, mapping, seg);
    compressSegmentation(seg);
end

