function job = aggloGridSeachSimple(p)
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
borderMeta = load([p.saveFolder 'globalBorder.mat'], 'borderSize', 'borderCoM');
toc;

% do grid search on these parameter values
borderSize = [25 50 100 200 300 500]; % border below this size are excluded from dendrite subgraph during inital agglomeration
segmentSize = [25 50 100 300 500 700 1000]; % segment below ...
probThreshold = [0.94 0.95 0.96 0.97 0.98 0.99]; % threshold on neurite continuity probability for CC
sizeThreshold = 100; % threshold on final agglomerate size in voxels (single segments not collected!)

Util.log('create job input parameters')
inputArguments = inputArgumentsFromParameterSets(borderSize, segmentSize, probThreshold, sizeThreshold);
inputArgumentFilename = arrayfun(@(x) fullfile(outputFolder, num2str(x,'%.5i'),'.mat'), 1:size(inputArguments,1),'uni',0)';
inputArguments = mat2cell(inputArguments, ones(size(inputArguments,1),1), size(inputArguments,2));
inputArguments = cellfun(@(x,y)[num2cell(x) y], inputArguments, inputArgumentFilename, 'uni', 0);

Util.log('Submitting jobs...')
job = Cluster.startJob( ...
    @doAgglomeration, inputArguments, ...
    'sharedInputs', {graph, segmentMeta, borderMeta, p, info}, ...
    'sharedInputsLocation', [6,7,8,9,10], ...
    'name', 'gridSearch', ...
    'cluster', {'memory', 12, ...
        'time', '24:00:00', 'priority', 100});
end

function [agglos, agglosSize, agglosEdges] = doAgglomeration(...
            borderSizeThr,segmentSizeThr, probThreshold, sizeThreshold, outputFile,...
            graph, segmentMeta, borderMeta, p, info)
    display(['Cut graph at border size:' num2str(borderSizeThr), 'segment size:' num2str(segmentSizeThr)]);
    graphCut = connectEM.cutGraphSimple(p, graph, segmentMeta, borderMeta, borderSizeThr, segmentSizeThr);
    display(['Performing agglomeration on graph with thr prob:' num2str(probThreshold), 'agglo size:' num2str(sizeThreshold)]);
    [agglos, agglosSize, agglosEdges] = connectEM.partitionSortAndKeepOnlyLarge(graphCut, segmentMeta, probThreshold, sizeThreshold);
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

