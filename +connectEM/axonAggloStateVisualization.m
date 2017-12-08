function axonAggloStateVisualization()
    % Given a superagglo axon state perform the following steps:
    % 1) Collect small axon superagglos along flight path of large superagglo
    % 2) Determine which of these augmented superagglos overlap with any axon GT
    % 3) Write gallery of GT skeletons and overlapping superagglos
    % 4) Write same to nmls for visualization in webKnossos
    
    %% Settings
    nodeEvidenceNeighbourhood = 3;
    minNodeEvidence = 54;
    minSegmentOverlap = 3;
    
    rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI/';
    axonFile = fullfile(rootDir, 'aggloState', 'axons_13_a.mat');
    
    thisDir = fileparts(mfilename('fullpath'));
    axonGtDir = fullfile(thisDir, 'evaluationData', 'new_axon_gt_ROI2017');
    cacheDir = '/tmpscratch/amotta/l4/2017-12-08-axons-13-visualization';
    
   [~, axonName, axonExt] = fileparts(axonFile);
    axonFlightFile = fullfile(cacheDir, strcat(axonName, axonExt));

    %% Load & modify parameter struct to use WKW segmentation
    param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
    param = param.p;
    param.seg = struct;
    param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
    param.seg.backend = 'wkwrap';

    if ~exist(cacheDir, 'dir')
        mkdir(cacheDir)
    end

    % Step 0.5: Load superagglos
    disp('Loading axon superagglos:');
    tic;
    axons = load(axonFile);
    toc;

    % (slightly modified from +connectEM.buildAxonAgglomerates.m)
    if exist(axonFlightFile, 'file')
        disp('Loading segment ids along flights paths from cache:');
        tic;
        load(axonFlightFile, 'axonFlights', 'axonFlightsMeta');
        toc;
    else
        disp('Looking up segment ids along flight paths:');
        tic;
        [axonFlights, axonFlightsMeta] = Superagglos.getFlightPathSegIds(param, axons.axons(axons.indBigAxons), nodeEvidenceNeighbourhood);
        Util.saveStruct(axonFlightFile, struct('axonFlights', axonFlights, 'axonFlightsMeta', axonFlightsMeta));
        toc;
    end
    
    disp('Picking up small superagglos along flight path of large superagglos:');
    tic;
    fullAxonAgglos = pickupSmallAxonAgglomerates(param, axons, axonFlights, minNodeEvidence);
    toc;

    disp('Reading GT skeletons:');
    tic;
    gtFiles = dir(fullfile(axonGtDir, '*.nml'));
    skel = arrayfun(@(x)skeleton(fullfile(axonGtDir, x.name)), gtFiles, 'uni', 0);
    toc;

    disp('Writing overlap PDF and skeletons');
    tic;
    segmentMeta = load(fullfile(param.saveFolder, 'segmentMeta.mat'), 'point');
    skelToAgglos = L4.Agglo.aggloSkelOverlap(skel, param, fullAxonAgglos);
    ovSkels = connectEM.skelOverlapGallery(cacheDir, skel, fullAxonAgglos, skelToAgglos, minSegmentOverlap, segmentMeta.point');
    for i = 1:length(ovSkels)
        ovSkels{i}.write(fullfile(cacheDir, ['gtSkel' num2str(i, '%.2i') '.nml']));
    end
    toc;

end

function axonAgglos = pickupSmallAxonAgglomerates(param, axons, axonFlights, minEvidence)
    % Add small (< 5 micron) axon agglomerates to large on if they have enough overlap
    % from +connectEM.buildAxonAgglomerates.m

    % small agglomerates to pick up
    axonAgglosSmall = axons.axons(~axons.indBigAxons);
    axonAgglosSmall = Superagglos.getSegIds(axonAgglosSmall);

    maxSegId = Seg.Global.getMaxSegId(param);
    axonAgglosLUT = Agglo.buildLUT(maxSegId, axonAgglosSmall);
    axonAgglosLUT = cat(1, 0, axonAgglosLUT(:));

    axonFlights.smallAxonId = axonAgglosLUT(1 + axonFlights.segId);
    axonFlights(~axonFlights.smallAxonId, :) = [];

    % pool evidence over agglomerates
    [axonOverlap, ~, axonEvidence] = unique(axonFlights(:, {'aggloId', 'flightId', 'smallAxonId'}), 'rows');
    axonOverlap.evidence = accumarray(axonEvidence, 1);

    % discard overlaps below evidence threshold
    axonOverlap = sortrows(axonOverlap, 'evidence', 'descend');
    axonOverlap(axonOverlap.evidence < minEvidence, :) = [];

    % assign to axon with largest evidence
    [~, uniRows] = unique(axonOverlap.smallAxonId, 'stable');
    axonOverlap = axonOverlap(uniRows, :);

    % build final agglomerates
    axonAgglosLarge = axons.axons(axons.indBigAxons);
    axonAgglosLarge = Superagglos.getSegIds(axonAgglosLarge);

    pickedUpAgglos = accumarray( ...
        axonOverlap.aggloId, axonOverlap.smallAxonId, size(axonAgglosLarge), ...
        @(r) {cat(1, axonAgglosSmall{r})}, {zeros(0, 1)});

    axonAgglos = cellfun( ...
        @vertcat, axonAgglosLarge, pickedUpAgglos, ...
        'UniformOutput', false);
end

