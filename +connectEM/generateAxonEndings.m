function generateAxonEndings(p, graph, segmentMeta, borderMeta, directionality, axons, outputFolder, options)
addpath(genpath('/u/scchr/Repositories/auxiliaryMethods/'));
addpath(genpath('/u/scchr/Repositories/pipeline/'));
addpath(genpath('/u/scchr/Repositories/nuclear_pores/'));

    options.latentScore = 0.5;
    options.segDirScore = 0.8;
    options.border = [3000; -3000];
    options.writeTasksToFile = true;
    options.boundingBoxForTasks = false;

    load('/gaba/scratch/mberning/axonQueryGeneration/beforeQueryGeneration.mat' , 'directionality')%, 'axonsNew')
    load('destinationOfModifiedAxonRepresentation')
    
    
%     axons = axonsNew;
    load('/tmpscratch/scchr/AxonEndings/Meta.mat','borderMeta')
    load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
    


    % Find all endings (given criteria on latent, directionality, orientation along z-axis)
    idxDirectional = cellfun(@(x)x(:,1) > options.latentScore, directionality.latent, 'uni', 0);
    idxEnding = cellfun(@(x)abs(x) > options.segDirScore, directionality.scores, 'uni', 0);
    idxAll = cellfun(@(x,y)find(x&y), idxDirectional, idxEnding, 'uni', 0);
    % Keep only largest score in each agglomerate for now
    nrCanidates = cellfun(@numel, idxAll);
    idxCanidateFound = nrCanidates > 0;
    AggloIDs=[1:length(idxAll)]';
    mapping=AggloIDs(idxCanidateFound);
    
    % clustering on left candidates
    thisBorderIdx = cellfun(@(x,y)x(y), directionality.borderIdx(idxCanidateFound), idxAll(idxCanidateFound),'uni',0);
    borderPositions = cellfun(@(x) borderMeta.borderCoM(x,:), thisBorderIdx,'uni',0);
    borderPositions_0 = cellfun(@(x) borderMeta.borderCoM(x,:), directionality.borderIdx(idxCanidateFound),'uni',0);

    clear borderMeta
    
    distance_cutoff = 600;
    T = cellfun(@(x) AxonEndings.clusterSynapticInterfaces(x,distance_cutoff,[11.24 11.24 28]), borderPositions, 'uni', 0);
      
    
    