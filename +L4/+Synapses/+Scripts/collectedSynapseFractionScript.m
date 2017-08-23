%% script to determine the already fraction of detected synapses
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%% load graph data

stats = struct(); % stores result
m = load('E:\workspace\data\backed\20170217_ROI\globalEdges.mat', 'edges');
edges = m.edges;
m = load('E:\workspace\data\backed\20170217_ROI\globalSynScores.mat');
synScores = nan(size(edges));
synScores(m.edgeIdx, :) = m.synScores;
synT = -1.67;
stats.synT = synT;
synEdgeIdx = max(synScores, [], 2) > -1.67;
m = load('E:\workspace\data\backed\20170217_ROI\graph.mat');
prob = m.prob(~isnan(m.borderIdx));

% discard synEdges with high continuity prob
probT = 0.5;
stats.probT = probT;
synEdgeIdx = synEdgeIdx & prob < 0.5;

% margin at border of dataset
bbox = [ 129, 5574; 129, 8509; 129, 3414];
margin = round(3000./[11.24; 11.24; 28]);
bboxM = bsxfun(@plus, bbox, [margin, -margin]);
m = load('E:\workspace\data\backed\20170217_ROI\globalBorder.mat', ...
    'borderCoM');
borderCom = m.borderCoM;
isAtBorder = ~Util.isInBBox(borderCom, bboxM);

stats.margin = margin; % nm

%% get axon agglos > 5um

m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
m.isAbove5um(m.toIgnoreIdx) = false;
axonsPostQuery = m.axonsPostQuery(m.isAbove5um);
axonIds = cell2mat(axonsPostQuery);

%% find fraction that overlap (>5um agglos)

% whole dataset
synEdges = edges(synEdgeIdx, :);
collectedSynEdges_5um = any(ismember(synEdges, axonIds), 2);
stats.fraction_5um = sum(collectedSynEdges_5um)./length(collectedSynEdges_5um);
stats.total_5um = sum(collectedSynEdges_5um);

% not at border
synEdges = edges(synEdgeIdx & ~isAtBorder, :);
collectedSynEdges_5um = any(ismember(synEdges, axonIds), 2);
stats.fraction_5um_margin = sum(collectedSynEdges_5um)./length(collectedSynEdges_5um);
stats.total_5um_margin = sum(collectedSynEdges_5um);

%% get all axon agglos

m = load('E:\workspace\data\backed\L4\FocusEM\postQueryAgglomerates.mat');
axonsPostQuery = m.axonsPostQuery;
axonIds = cell2mat(axonsPostQuery);

%% find fraction that overlap

synEdges = edges(synEdgeIdx, :);
collectedSynEdges = any(ismember(synEdges, axonIds), 2);
stats.fraction = sum(collectedSynEdges)./length(collectedSynEdges);
stats.total = sum(collectedSynEdges);

% not at border
synEdges = edges(synEdgeIdx & ~isAtBorder, :);
collectedSynEdges_5um = any(ismember(synEdges, axonIds), 2);
stats.fraction_margin = sum(collectedSynEdges_5um)./length(collectedSynEdges_5um);
stats.total_margin = sum(collectedSynEdges_5um);

%% display results
disp(stats)
% save('E:\workspace\data\backed\L4\Synapses\collectedSynapseFraction.mat', 'stats')