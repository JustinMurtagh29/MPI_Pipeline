function generateAxonEndingInputData(param)
    % generateAxonEndingInputData(param)
    %   Generates a `endingInputData.mat` file in the pipeline directory
    %   which contains all the data necessary for the ending detection.
    %
    % Written by
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>

    % Parameters
    minAxonLength = 5000; % in nm
    bboxDist = 1000; % in nm

    % Directory with input / output data
    dataDir = fullfile(param.saveFolder, 'aggloState');
    intermediateFile = fullfile(dataDir, 'axonDirectionality.mat');
    outFile = fullfile(dataDir, 'axonEndingInputData.mat');

    % Load data
   [graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);

    % Load state of axon agglomeration
    axons = load(fullfile(dataDir, 'axons_04.mat'));
    axons = axons.axons;

    % Convert to old-school agglomerates
    axonAgglos = arrayfun( ...
        @Agglo.fromSuperAgglo, axons, 'UniformOutput', false);

    % Calculate axon directionality
    options = struct;
    options.voxelSize = param.raw.voxelSize;
    options.bboxDist = bboxDist;

    directionality = connectEM.calculateDirectionalityOfAgglomerates( ...
        axonAgglos, graph, segmentMeta, borderMeta, globalSegmentPCA, options);
    save(intermediateFile, 'directionality', '-v7.3');

    % Convert borders to nm
    borderMeta.borderCoM = bsxfun( ...
        @times, double(borderMeta.borderCoM), param.raw.voxelSize);

    axonMaskFunc = @(b) Util.isMaxPdistAbove( ...
        borderMeta.borderCoM(b, :), minAxonLength);
    axonIds = find(cellfun(axonMaskFunc, directionality.borderIdx));

    % Save results
    out = struct;
    out.axons = axons(axonIds);
    out.axonIds = axonIds;
    out.directionality = structfun( ...
        @(x) x(axonIds), directionality, 'UniformOutput', false);
    out.gitInfo = Util.gitInfo();

    Util.saveStruct(outFile, out);

end

