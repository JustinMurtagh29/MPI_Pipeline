function axonAggloStateVisualization()
    % Given a superagglo axon state perform the following steps:
    % 1) Collect small axon superagglos along flight path of large superagglo
    % 2) Determine which of these augmented superagglos overlap with any axon GT
    % 3) Write gallery of GT skeletons and overlapping superagglos
    % 4) Write same to nmls for visualization in webKnossos
    
    % HACK(amotta):
    %   This function depends on code from Benedikt's repository in order
    %   to generate the overlap skeletons and gallery. His repository used
    %   to be a submodule of the pipeline. Since this is no longer true, we
    %   temporarily add his repository to the search path.
    oldPath = addpath('/gaba/u/amotta/code/benedikt');
    cleanup = onCleanup(@() path(oldPath));
    
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
   [completedAgglos, pickUpIds] = ...
        pickupSmallAxonAgglomerates( ...
            param, axons, axonFlights, minNodeEvidence);
    toc;

    disp('Reading GT skeletons:');
    tic;
    gtFiles = dir(fullfile(axonGtDir, '*.nml'));
    skel = arrayfun( ...
        @(x) skeleton(fullfile(axonGtDir, x.name)), ...
        gtFiles, 'UniformOutput', false);
    toc;

    disp('Writing overlap PDF and skeletons');
    tic;
    
    % determine overlap between skeleton tracings and completed agglos
    skelToAgglos = L4.Agglo.aggloSkelOverlap(skel, param, completedAgglos);
    
    % replace completed agglo ID with the IDs of all its components
    skelToAgglos = cellfun( ...
        @(o) horzcat( ...
            cat(1, pickUpIds{o(:, 1)}), ...
            repelem(o(:, 2), cellfun( ...
                @numel, pickUpIds(o(:, 1))), 1)), ...
        skelToAgglos, 'UniformOutput', false);
    
    ovSkels = connectEM.skelOverlapGallery( ...
        cacheDir, skel, axons.axons, skelToAgglos, minSegmentOverlap);
    
    for i = 1:length(ovSkels)
        ovSkels{i}.parameters.experiment.name = param.experimentName;
        ovSkels{i}.write(fullfile(cacheDir, ['gtSkel' num2str(i, '%.2i') '.nml']));
    end
    toc;
end

function [completedAgglos, pickUpIds] = ...
        pickupSmallAxonAgglomerates(param, axons, axonFlights, minEvidence)
    % Add small (< 5 micron) axon agglomerates to large on if they have enough overlap
    % from +connectEM.buildAxonAgglomerates.m
    
    largeIds = find(axons.indBigAxons(:));
    smallIds = find(~axons.indBigAxons(:));
    agglos = Superagglos.getSegIds(axons.axons);
    
    %% determine overlap
    maxSegId = Seg.Global.getMaxSegId(param);
    smallLUT = Agglo.buildLUT(maxSegId, agglos(smallIds));
    smallLUT = cat(1, 0, smallLUT(:));

    axonFlights.smallId = smallLUT(1 + axonFlights.segId);
    axonFlights(~axonFlights.smallId, :) = [];
    axonFlights.smallId = smallIds(axonFlights.smallId);

    % pool evidence over agglomerates
   [axonOverlap, ~, axonEvidence] = unique( ...
       axonFlights(:, {'aggloId', 'flightId', 'smallId'}), 'rows');
    axonOverlap.evidence = accumarray(axonEvidence, 1);

    % discard overlaps below evidence threshold
    axonOverlap = sortrows(axonOverlap, 'evidence', 'descend');
    axonOverlap(axonOverlap.evidence < minEvidence, :) = [];

    % assign to axon with largest evidence
   [~, uniRows] = unique(axonOverlap.smallId, 'stable');
    axonOverlap = axonOverlap(uniRows, :);

    %% build final agglomerates
    pickUpIds = accumarray( ...
        axonOverlap.aggloId, axonOverlap.smallId, ...
       	size(largeIds), @(ids) {ids}, {zeros(0, 1)});
    
    pickUpIds = cellfun( ...
        @vertcat, num2cell(largeIds), pickUpIds, 'UniformOutput', false);
    completedAgglos = cellfun( ...
        @(ids) cat(1, agglos{ids}), pickUpIds, 'UniformOutput', false);
end
