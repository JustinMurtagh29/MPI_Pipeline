connName = 'connectome_ax18a_deWC01wSp';
mainOutDir = ['/tmpscratch/mbeining/isos/',connName];
isoParamAx = {'smoothWidth',5,'smoothSizeHalf',5,'reduce',0.1};  % heiko suggests 9 , 9 , 0.2
isoParamMy = {'smoothWidth',2,'smoothSizeHalf',2,'reduce',0.15};

%% load stuff
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');
connDir = fullfile([p.saveFolder], 'connectomeState');

if ~exist('graph','var') || ~all(isfield(graph,{'edges','prob','borderIdx'}))
    graph = load(fullfile(p.saveFolder, 'graphNew.mat'),'edges','prob','borderIdx');
end
if  ~all(isfield(graph,{'neighbours','neighBorderIdx'}))
    [graph.neighbours, neighboursIdx] = Graph.edges2Neighbors(graph.edges);
    graph.neighBorderIdx = cellfun(@(x)graph.borderIdx(x), neighboursIdx, 'uni', 0);
    clear neighboursIdx
end
if ~exist('borderMeta','var')
    borderMeta = load(fullfile(p.saveFolder, 'globalBorder.mat'), 'borderSize', 'borderCoM');
end
if ~exist('heuristics','var')
    heuristics = load(fullfile(p.saveFolder, 'heuristicResult.mat'),'myelinScore');
end
connFile = fullfile(connDir, sprintf('%s.mat', connName));

% loading connectome data
conn = load(connFile);

interSynFile = sprintf('%s_intersynapse.mat', connName);
interSynFile = fullfile([connDir], interSynFile);
interSyn = load(interSynFile);

disp('all data loaded')
%%
% get classification thalamocortical and myelinated
isThalamocortical = connectEM.Axon.detectThalamocorticals(conn, interSyn);
[ myelinFracAx,~,~,myelinNeighSegments] = connectEM.calculateSurfaceMyelinScore( conn.axons, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class

disp('Done calculating TC/myelinated axons, creating isos')
%% export TC axons including myelin iso
Visualization.exportAggloToAmira(p, conn.axons(isThalamocortical), fullfile(mainOutDir,'TC'), isoParamAx{:})
Visualization.exportAggloToAmira(p, myelinNeighSegments(isThalamocortical), fullfile(mainOutDir,'TCmyelin'), isoParamMy{:})


%% export myelinated axons including myelin iso, excluding tc axons
myAxNoTC = ~isThalamocortical & myelinFacAx > 0.2;
Visualization.exportAggloToAmira(p, conn.axons(myAxNoTC), fullfile(mainOutDir,'AxMy'), isoParamAx{:})
Visualization.exportAggloToAmira(p, myelinNeighSegments(myAxNoTC), fullfile(mainOutDir,'AxMymyelin'), isoParamMy{:})
