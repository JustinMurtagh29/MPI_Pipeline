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

    % Load data
   [graph, segmentMeta, borderMeta, globalSegmentPCA] = ...
       connectEM.loadAllSegmentationData(param);

    % Load state of axon agglomeration
%     axonsNew = '/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/metricsFinal.mat';
%     axonsNew = load(axonsNew, 'axonsNew');
%     axonsNew = axonsNew.axonsNew;
    m = load(fullfile(param.saveFolder, 'aggloState/', 'axons_04.mat'));
    axonsNew = m.axons;
    
    % Calculate axon directionality
    options = struct;
    options.voxelSize = param.raw.voxelSize;
    options.bboxDist = bboxDist;
    
    directionality = connectEM.calculateDirectionalityOfAgglomerates( ...
        axonsNew, graph, segmentMeta, borderMeta, globalSegmentPCA, options);

    % Calculate maximum border-to-border size within agglomerates
    calculateLength = @(borderIds) max(pdist(bsxfun( ...
        @times, double(borderMeta.borderCoM(borderIds, :)), param.raw.voxelSize)));
    axonLength = cellfun(calculateLength, directionality.borderIdx, 'uni', 0);

    % Keep only long enough axons
    idxKeep = ~cellfun(@isempty, axonLength);
    idxKeep(idxKeep) = cellfun(@(x) x > minAxonLength, axonLength(idxKeep));
    
    % Save results
    out = struct;
    out.idxKeep = idxKeep;
    out.axonsNew = axonsNew(idxKeep);
    out.directionality = structfun( ...
        @(x) x(idxKeep), directionality, 'UniformOutput', false);
    out.axonLength = cell2mat(axonLength(idxKeep));
    out.gitInfo = Util.gitInfo();
    
    outFile = fullfile(param.saveFolder, 'aggloState/', 'axonEndingInputData.mat');
    Util.saveStruct(outFile, out);
end
