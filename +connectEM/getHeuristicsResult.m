function getHeuritisticResult(p)
    % Collect results of vessel, nuclei and mylein detection

        % Get vessel scores
        binaryMap(1).root = [p.saveFolder,'NucleiAndVesselDetection/segmentationVessel/1/';
        binaryMap(1).prefix = 'vessel_mag1'; 
        binaryMap(1).segId = 1;
        % Get myelin scores
        binaryMap(2).root = '/gaba/u/mberning/wkCubes/2012-09-28_ex145_07x2_ROI2016_vessel/segmentation/1/';
        binaryMap(2).prefix = '2012-09-28_ex145_07x2_ROI2016_vessel_mag1'; 
        binaryMap(2).segId = 3;
        % Get nuclei scores
        binaryMap(3).root = [p.saveFolder,'NucleiAndVesselDetection/segmentationNuclei/1/'];
        binaryMap(3).prefix = 'nuclei_mag1'; 
        binaryMap(3).segId = 2;


    for cube=1:numel(p.local)
        bbox = p.local(cube).bboxSmall;
        inputCell{cube} = {bbox};
    end
    functionH = @connectEM.getSegmentHeuristicsScore;

    % Execute and retrieve results
    seg = p.seg;
    cluster = Cluster.getCluster('-l h_vmem=12G');
    job = Cluster.startJob(functionH, inputCell, ...
        'name', 'heuristicLookup', 'sharedInputs', {seg binaryMap}, ...
        'sharedInputsLocation', [1 2], 'cluster', cluster, ...
        'numOutputs', 2);
    Cluster.waitForJob(job);
    scores = fetchOutputs(job);

    % Reformatting for save file
    segIds = cat(1, scores{:,2});
    sc = cat(1, scores{:,1}); 
    vesselScore = cat(1, sc{:,1});
    myelinScore = cat(1, sc{:,2});
    nucleiScore = cat(1, sc{:,3});

    Util.save([p.saveFolder 'heuristicResult.mat'], segIds, vesselScore, myelinScore, nucleiScore);
end

function scoresSub = extractScores(scores, idx)
    temp = cellfun(@(x)x{idx}, scores(:,1), 'uni', 0);
    scoresSub = cell2mat(temp);
end
