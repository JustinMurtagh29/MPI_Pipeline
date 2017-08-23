function endingDetection(p, borderMeta, directionality, axonsNew)

    options.latentScore = 0.5;
    options.segDirScore = 0.8;
    options.border = [3000; -3000];
    options.writeTasksToFile = true;
    options.boundingBoxForTasks = false;

    for i=1:length(axonsNew)
        axons{i,1} = axonsNew(i).nodes(:,4);
    end
    axons = cellfun(@(x)x(~isnan(x)),axons,'uni',0);
    
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    
    % Find all endings (given criteria on latent, directionality, orientation along z-axis)
    idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, directionality.latent, 'uni', 0);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, directionality.scores, 'uni', 0);
    idxAll = cellfun(@(x,y)find(x&y), idxDirectional, idxEnding, 'uni', 0);
    % Keep only largest score in each agglomerate for now
    nrCanidates = cellfun(@numel, idxAll);
    idxCanidateFound = nrCanidates > 0;
    
    % clustering on left candidates
    thisBorderIdx = cellfun(@(x,y)x(y), directionality.borderIdx(idxCanidateFound), idxAll(idxCanidateFound),'uni',0);
    % borderPositions of all candidates fullfilling the foregoing
    % conditions
    borderPositions = cellfun(@(x) borderMeta.borderCoM(x,:), thisBorderIdx,'uni',0);
    % borderPositions of the whole agglomerate (except those agglos without
    % a single candidate for an ending)
    borderPositions_0 = cellfun(@(x) borderMeta.borderCoM(x,:), directionality.borderIdx(idxCanidateFound),'uni',0);

    clear borderMeta
    
    distance_cutoff = 600;
    T = cellfun(@(x) AxonEndings.clusterSynapticInterfaces(x,distance_cutoff,[11.24 11.24 28]), borderPositions, 'uni', 0);
    
    save([p.savefolder 'AxonEndingDetectionResults.mat'])
      
    
    