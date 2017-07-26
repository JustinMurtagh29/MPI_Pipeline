% Did not have access to /tmpscratch/kboerg/aggloGridSearch6/6_01_00046/, so I redid with some more code moved here for completeness

%{
% Find arguments used in directionality based agglomeration, see aggloGridSearchDir
b = load('/gaba/scratch/mberning/aggloGridSearch6/parameters.mat');
options.latentScore = b.inputArgumentsAxons(46,1);
options.segDirScore = b.inputArgumentsAxons(46,2);
options.neuriCScore = b.inputArgumentsAxons(46,3);
options.borderSize = b.inputArgumentsAxons(46,4);
options.axonScore = b.inputArgumentsAxons(46,5);
options.sourceSize = 2000;
options.recursionSteps = 10;
options.minSize = 100;
options.bboxDist = 1000;
options.voxelSize = [11.24 11.24 28];
inputArguments = {options '/tmpscratch/mberning/axonQueryResults/dirGridSearchRedo2/'};
connectEM.axonDirectionalityBasedGrowing(inputArguments{:});
clear b options inputArguments;

% Make sure the directionality based result is the same as before
a = load('/gaba/scratch/mberning/aggloGridSearch6/6_01_00046/10.mat');
b = load('/tmpscratch/mberning/axonQueryResults/dirGridSearchRedo2/10.mat');
assert(all(cellfun(@(x,y)all(x==y), a.axonsNew, b.axonsNew)));
clear a b;
%}

% Load graph
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
[graph, segmentMeta, borderMeta, globalSegmentPCA] = connectEM.loadAllSegmentationData(p);

% Get edges switched on during directionality based growing
bordersIdx = [];
for idx = 1:10
    temp = load(['/tmpscratch/mberning/axonQueryResults/dirGridSearchRedo2/' sprintf('%0.2u',idx)]);
    bordersIdx = [bordersIdx; cell2mat(temp.edgesToStore)];
end
edgesIdx = find(ismember(graph.borderIdx, bordersIdx));
% Add edges from inital CC based agglomeration also used as staring point for directionality based agglomeration
startAgglo = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat');
edgesGTall = [graph.edges(edgesIdx,:); startAgglo.axonEdges];
% Also add edges added during garbage collection
[axonsNew, ~, garbageEdges] = connectEM.garbageCollection(graph, segmentMeta, startAgglo.axons, startAgglo.dendrites, startAgglo.heuristics.mapping)
assert(all(cellfun(@(x,y)all(x==y), axonsNew, startAgglo.axonsFinal)));
edgesGTall = [graph.edges(edgesIdx,:); startAgglo.axonEdges; garbageEdges];
clear bordersIdx idx temp edgesIdx startAgglo axonsNew garbageEdges;
save('/tmpscratch/mberning/edgesGTall.mat', 'edgesGTall');

% Load eqClass based on _a and _b axon wK query projects
load('/tmpscratch/mberning/axonQueryResults/postQueryState.mat', ...
    'eqClassCCfull', 'ff', 'axons', 'startAgglo', 'endAgglo', 'idxGood');

% Generate a representation of links all queries make
edgesBetweenAgglos = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
queryLinks = cat(2, repelem(1:numel(edgesBetweenAgglos), cellfun(@(x)size(x,1), edgesBetweenAgglos))', sort(cat(1, edgesBetweenAgglos{:}),2));
% Only keep unique querry links
[~,idxA] = unique(queryLinks(:,2:3));
queryLinks = queryLinks(idxA,:);
clear idxA;

for idx_eqClass = 1000:1100
    tic;
    currentEC = eqClassCCfull{idx_eqClass};
    currentAgglo = axons(currentEC);
    assert(all(ismember(cat(1,currentAgglo{:}), edgesGTall)));
    % Get all segment positions in this eqClass and all edges between them
    aggloNodes = cellfun(@(x)segmentMeta.point(:, x)', currentAgglo, 'uni', 0);
    [locA, locB] = cellfun(@(x)ismember(edgesGTall, x), currentAgglo, 'uni', 0);
    aggloEdges = cellfun(@(x,y)y(all(x, 2), :), locA, locB, 'uni', 0);
    % Get all query node position and respective edges
    queryLinkIdx = queryLinks(all(ismember(queryLinks(:,2:3), currentEC),2),1);
    queriesNodes = [];
    queriesEdges = [];
    for idxQuery=1:length(queryLinkIdx)
        skel = skeleton(ff.filenames{queryLinkIdx(idxQuery)});
        queriesNodes{idxQuery,1} = skel.nodes{1}(:,1:3);
        queriesEdges{idxQuery,1} = skel.edges{1};
    end
    % Globalize representation of nodes and edges
    aggloNrNodes = cellfun(@(x)size(x,1), aggloNodes);
    queriesNrNodes = cellfun(@(x)size(x,1), queriesNodes);
    nodeOffset = cumsum(cat(1, 0, aggloNrNodes, queriesNrNodes(1:end-1)));
    nodes = cat(1, aggloNodes{:}, queriesNodes{:});
    edges = cellfun(@plus, cat(1, aggloEdges, queriesEdges), num2cell(nodeOffset), 'uni', 0);
    edges = cat(1, edges{:});
    % Join query and agglo based representation by single edge at end of query to closest node in connecting agglo
    for idxQuery=1:length(queryLinkIdx)
        idxDegreeBelow2 = find(histc(queriesEdges{idxQuery}(:), unique(queriesEdges{idxQuery})) < 2);
        assert(numel(idxDegreeBelow2) == 2);
        for i=1:numel(idxDegreeBelow2)
            % Left TODO: Find good way of making connections
            queryLinks(queryLinks(:,1) == queryLinkIdx)
        end
    end
    % Detect chiasmata
    outputFolder = ['/tmpscratch/mberning/axonQueryResults/chiasmataSkeletons/' num2str(idx_agglo, '%.4i') '/'];
    detectChiasmata(nodes, edges, true, outputFolder);
    time = toc;
    display(['Skeleton finished, nodes: ' num2str(size(nodes,1)) ', time: ' num2str(time) ', time/node: ' num2str(size(nodes,1)./time)]);
end

