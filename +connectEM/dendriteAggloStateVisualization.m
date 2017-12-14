function dendriteAggloStateVisualization()
    % Based on `axonAggloStateVisualization`
    %
    % Modified by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    
    % HACK(amotta):
    %   This function depends on code from Benedikt's repository in order
    %   to generate the overlap skeletons and gallery. His repository used
    %   to be a submodule of the pipeline. Since this is no longer true, we
    %   temporarily add his repository to the search path.
    oldPath = addpath('/gaba/u/amotta/code/benedikt');
    cleanup = onCleanup(@() path(oldPath));
    
    %% Settings
    minSegmentOverlap = 3;
    
    rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI/';
    dendFile = fullfile(rootDir, 'aggloState', 'dendrites_16.mat');
    
    dendGtDir = fullfile( ...
        fileparts(mfilename('fullpath')), ...
        'evaluationData', 'dendrite_gt_spine_head_seeded');
    outputDir = '/tmpscratch/amotta/l4/2017-12-14-dendrites-16-visualization';
    assert(exist(outputDir, 'dir') ~= 0);
    
    %% Load & modify parameter struct to use WKW segmentation
    param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
    param = param.p;
    param.seg = struct;
    param.seg.root = '/tmpscratch/amotta/l4/2012-09-28_ex145_07x2_ROI2017/segmentation/1';
    param.seg.backend = 'wkwrap';

    % Step 0.5: Load superagglos
    disp('Loading superagglos');
    dends = load(dendFile);

    disp('Reading GT skeletons');
    gtFiles = dir(fullfile(dendGtDir, '*.nml'));
    skel = arrayfun( ...
        @(x) skeleton(fullfile(dendGtDir, x.name)), ...
        gtFiles, 'UniformOutput', false);

    % determine overlap between skeleton tracings and completed agglos
    disp('Determining skeleton-agglomerate overlaps');
    dendAgglos = Superagglos.getSegIds(dends.dendrites);
    skelToAgglos = L4.Agglo.aggloSkelOverlap(skel, param, dendAgglos);
    
    disp('Writing overlap PDFs');
    ovSkels = connectEM.skelOverlapGallery( ...
        outputDir, skel, dends.dendrites, skelToAgglos, minSegmentOverlap);
    
    disp('Writing overlap NMLs');
    for i = 1:length(ovSkels)
        ovSkels{i}.parameters.experiment.name = param.experimentName;
        ovSkels{i}.write(fullfile(outputDir, ['gtSkel' num2str(i, '%.2i') '.nml']));
    end
end
