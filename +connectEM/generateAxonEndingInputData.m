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

    % Calculate maximum border-to-border size within agglomerates
    calculateLength = @(borderIds) max(pdist(bsxfun( ...
        @times, double(borderMeta.borderCoM(borderIds, :)), param.raw.voxelSize)));
    axonLengths = cellfun(calculateLength, directionality.borderIdx, 'uni', 0);

    % Keep only long enough axons
    indBigAxons = ~cellfun(@isempty, axonLengths);
    indBigAxons(indBigAxons) = cellfun( ...
        @(x) x > minAxonLength, axonLengths(indBigAxons));
    
    % Save results
    out = struct;
    out.indBigAxons = indBigAxons;
    out.bigAxons = axons(indBigAxons);
    out.directionality = structfun( ...
        @(x) x(indBigAxons), directionality, 'UniformOutput', false);
    out.axonLengths = cell2mat(axonLengths(indBigAxons));
    out.gitInfo = Util.gitInfo();
    
    outFile = fullfile(dataDir, 'axonEndingInputData.mat');
    Util.saveStruct(outFile, out);
end
