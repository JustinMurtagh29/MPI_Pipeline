function getHeuristicsResult(p)
    % Collect results of vessel, nuclei and mylein detection

    % Get vessel scores
    binaryMap(1).root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
    binaryMap(1).prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
    binaryMap(1).segId = 1;
    % Get myelin scores
    binaryMap(2).root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
    binaryMap(2).prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1';
    binaryMap(2).segId = 3;
    % Get nuclei scores
    binaryMap(3).root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_nuclei/segmentation/1/';
    binaryMap(3).prefix = '2012-09-28_ex145_07x2_ROI2016_nuclei_mag1';
    binaryMap(3).segId = 2;
    
    info = Util.runInfo();
    info.param = rmfield(info.param, 'p');

    inputCell = arrayfun( ...
        @(l) {l.bboxSmall}, p.local, 'UniformOutput', false);
    functionH = @connectEM.getSegmentHeuristicsScore;

    % Execute and retrieve results
    job = Cluster.startJob(functionH, inputCell, ...
        'name', 'heuristicLookup', 'sharedInputs', {p.seg, binaryMap}, ...
        'cluster', {'memory', 12}, 'numOutputs', 2);
    Cluster.waitForJob(job);
    scores = fetchOutputs(job);

    % Reformatting for save file
    sc = cat(1, scores{:, 1}); 
    segIds = cat(1, scores{:, 2});
    vesselScore = cat(1, sc{:, 1});
    myelinScore = cat(1, sc{:, 2});
    nucleiScore = cat(1, sc{:, 3});
    
    saveFile = fullfile(p.saveFolder, 'heuristicResult.mat');
    
    Util.save( ...
        saveFile, info, segIds, ...
        vesselScore, myelinScore, nucleiScore);
end
