function [job, inputArguments] = aggloGridSeachSimple(p)
% use Manuel's grid search routine with simpler version
% not using heuristics, typeEM information
info = Util.runInfo();

outputFolder = [p.saveFolder datestr(clock,30) 'gridSearchAgglomeration/'];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

display('loading dependent variables...');
tic;
m = load([p.saveFolder 'allParameterWithSynapses.mat']);
p = m.p;
graph = load([p.saveFolder 'graph.mat']);
segmentMeta = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount', 'point', 'maxSegId');
segmentMeta.point = segmentMeta.point';
segmentMeta = connectEM.addSegmentClasInformation(segmentMeta);
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
toc;

% do grid search on these parameter values
borderSize = [25 50 100 200 300 500]; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSize = [25 50 100 300 500 700 1000]; % segment below ...
sizeThreshold = 100; % threshold on final agglomerate size in voxels (single segments not collected!)
dendriteProbThreshold = [0.5 0.6 0.7 0.8 0.9]; % restrictions of graph based on segment class proability
axonProbThreshold = [0.3 0.4 0.5 0.6 0.7 0.8 0.9]; % restrictions of graph based on segment class proability
probThresholdDendrite = [0.94 0.95 0.96 0.97 0.98 0.99]; % threshold on neurite continuity probability for CC
probThresholdAxon = [0.93 0.94 0.95 0.96 0.97 0.98]; % threshold on neurite continuity probability for CC

Util.log('create job input parameters')
inputArguments = inputArgumentsFromParameterSets(borderSize, segmentSize, sizeThreshold,...
                 dendriteProbThreshold, axonProbThreshold, probThresholdDendrite, probThresholdAxon);
inputArgumentFilename = arrayfun(@(x) fullfile(outputFolder, [num2str(x,'%.5i'),'.mat']), 1:size(inputArguments,1),'uni',0)';
inputArguments = mat2cell(inputArguments, ones(size(inputArguments,1),1), size(inputArguments,2));
inputArguments = cellfun(@(x,y)[num2cell(x) y], inputArguments, inputArgumentFilename, 'uni', 0);

Util.log('Submitting jobs...')
job = Cluster.startJob( ...
    @doAgglomeration, inputArguments, ...
    'sharedInputs', {graph, segmentMeta, borderMeta, p, info}, ...
    'sharedInputsLocation', [9,10,11,12,13], ...
    'name', 'gridSearch', ...
    'cluster', {'memory', 48, ... % 30 for human, 48 for Macaque
        'time', '24:00:00', 'priority', 100});
end

function [axons, dendrites] = doAgglomeration(...
            borderSizeThr,segmentSizeThr, sizeThreshold, dendClassProb, axonClassProb,...
            probThresholdDendrite, probThresholdAxon, outputFile, ...
            graph, segmentMeta, borderMeta, p, info)
    display(['Cut graph at border size:' num2str(borderSizeThr), 'segment size:' num2str(segmentSizeThr)]);
    graphCut = connectEM.cutGraphSimple(p, graph, segmentMeta, borderMeta, borderSizeThr, segmentSizeThr);
    Util.log('add class information to restrict graph')
    segmentMeta.isDendrite = segmentMeta.dendriteProb >= dendClassProb & ...
                                segmentMeta.gliaProb <= 0.85;
    segmentMeta.isAxon = segmentMeta.axonProb >= axonClassProb & ...
                                segmentMeta.gliaProb <= 0.85;

    Util.log('Generating subgraphs for axon and dendrite agglomeration:');
    idx = all(ismember(graphCut.edges, find(segmentMeta.isDendrite)), 2);
    graphCutDendrites.edges = graphCut.edges(idx,:);
    graphCutDendrites.prob = graphCut.prob(idx);
    [dendrites, dendriteSize, dendriteEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutDendrites,...
                                    segmentMeta, probThresholdDendrite, sizeThreshold);
    idx = all(ismember(graphCut.edges, find(segmentMeta.isAxon)), 2);
    graphCutAxons.edges = graphCut.edges(idx,:);
    graphCutAxons.prob = graphCut.prob(idx);
    [axons, axonSize, axonEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCutAxons,...
                                    segmentMeta, probThresholdAxon, sizeThreshold);
    display('Saving:');
    tic;
    % Lets save some more so that we can always understand whats happening
    clearvars borderMeta segmentMeta graph graphCut
    save(outputFile);
    toc;
end

function inputArguments = inputArgumentsFromParameterSets(a, b, c, d)

    [aG, bG, cG, dG ] = ndgrid(a, b, c, d);
    inputArguments = cat(2, aG(:), bG(:), cG(:), dG(:)); 

end

