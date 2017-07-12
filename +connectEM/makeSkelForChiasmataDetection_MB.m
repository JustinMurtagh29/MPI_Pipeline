% Did not have access to /tmpscratch/kboerg/aggloGridSearch6/6_01_00046/, so I redid with some more code moved here for completeness

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

% Get edges switched on during directionality based growing
bordersIdx = [];
for idx = 1:10
    temp = load(['/tmpscratch/mberning/axonQueryResults/dirGridSearchRedo/' sprintf('%0.2u',idx)]);
    bordersIdx = [bordersGT; cell2mat(temp.edgesToStore)];
end
edgesIdx = find(ismember(graph.borderIdx, bordersIdx));
% Add edges from inital CC based agglomeration also used as staring point for directionality based agglomeration
startAgglo = load('/gaba/scratch/mberning/aggloGridSearch/search05_00564.mat');
edgesGTall = [graph.edges(edgesIdx,:); startAgglo.axonEdges];
clear bordersIdx idx temp edgesIdx startAgglo;

% Load eqClass based on _a and _b axon wK query projects
load('/tmpscratch/mberning/axonQueryResults/postQueryState.mat', ...
    'eqClassCCfull', 'ff', 'axons', 'startAgglo', 'endAgglo', 'idxGood');

% Generate a representation of links all queries make
edgesBetweenAgglos = cellfun(@(x,y)combnk([x y], 2), startAgglo(idxGood), endAgglo(idxGood), 'uni', 0);
queryLinks = cat(2, repelem(1:numel(edgesBetweenAgglos), cellfun(@(x)size(x,1), edgesBetweenAgglos)), cat(1, edgesBetweenAgglos{:}));


for idx_eqClass = 1000:1100
    tic;
    currentEC = eqClassCCfull{idx_eqClass};
    currentAgglos = axons{currentEC};
    % Get all segment positions in this eqClass and all edges between them
    aggloNodes = cellfun(@(x)segmentMeta.point(:, x), currentAgglos, 'uni', 0);
    [locA, locB] = cellfun(@(x)ismember(edgesGTall, x), currentAgglos, 'uni', 0);
    aggloEdges = cellfun(@(x,y)y(all(x, 2), :), loxA, locB, 'uni', 0);
    % Get all query node position and respective edges
    queriesIdx = cellfun(@(x,y)any(ismember(x,currentEC))&any(ismember(y,currentEC)), startAgglo, endAgglo);

    queriesNodes = ff.nodes(idx);
    queriesEdges = cell(numel(queriesNodes),1);
    for idxQuery = 1:length(queriesNodes)
        queriesEdges{i} = cat(2, nrNodes+1:nrNodes+size(queriesNodes{i},1)-1, nrNodes+2:nrNodes+size(queriesNodes{i},1));

        nrNodes = nrNodes + size(queriesNodes{i},1);
    end
    % Join query and agglo based representation
    nrNodesOffsetA = cellfun(@(x)size(x,1)-size(aggloNodes{1},1), aggloNodes);
    nrNodesOffsetQ = cellfun(@(x)size(x,1)-size(queriesNodes{1},1), queriesNodes);
    aggloEdges = cellfun(@(x,y)x-y, aggloEdges, nrNodesOffsetA);
    queriesEdges = cellfun(@(x,y)x-y, queriesEdges, nrNodesOffsetQ);
    nodes = cat(1, aggloNodes{:}, queriesNodes{:});
    edges = cat(1, aggloEdges{:}, queriesEdges{:});

    % Detect chiasmata
    outputFolder = ['/tmpscratch/mberning/axonQueryResults/chiasmataSkeletons/' num2str(idx_agglo, '%.4i') '/'];
    detectChiasmata(nodes, edges, true, outputFolder);
    time = toc;
    display(['Skeleton finished, nodes: ' num2str(size(nodes,1)) ', time: ' num2str(time) ', time/node: ' num2str(size(nodes,1)./time)]);
end

