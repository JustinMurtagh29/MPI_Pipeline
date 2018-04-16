connFileName = 'connectome_ax18a_deWC01wSp_v4.mat';
mainOutDir = ['/tmpscratch/mbeining/isos/',connName];
isoParamAx = {'smoothWidth',5,'smoothSizeHalf',5,'reduce',0.1};  % heiko suggests 9 , 9 , 0.2
isoParamMy = {'smoothWidth',2,'smoothSizeHalf',2,'reduce',0.15};

%% load stuff
load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat');
p.backend = 'wkwrap';
disp('Parameters loaded');
outputFolder = fullfile(p.saveFolder, 'aggloState');

if ~exist('graph','var') || ~all(isfield(graph,{'edges','borderIdx'}))
    graph = load(fullfile(p.saveFolder, 'graphNew.mat'),'edges','borderIdx');
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

conn = connectEM.Connectome.load(p,fullfile(p.saveFolder,'connectomeState',connName));
load(fullfile(p.saveFolder,'aggloState','axons_18_a.mat'))

% map axons from connectome state to the ones from axon state in order to
% be able to apply the thalamocortical ids
[axonLUT,axonSegIds] = Superagglos.buildLUT(axons);
connAxonSegIds = cell2mat(conn.axons);
connAxonLUT = Agglo.buildLUT(numel(axonLUT),conn.axons);
[ismem,ind] = ismember(connAxonSegIds,axonSegIds);
connAxonIds = (accumarray(connAxonLUT(connAxonSegIds(ismem))',axonLUT(axonSegIds(ind(ismem)))',[],@mode));
axons = axons(connAxonIds);

disp('all data loaded')
%%
% get classification thalamocortical and myelinated
isThalamocortical = conn.axonMeta.isThalamocortical;
[ myelinFracAx,~,~,myelinNeighSegments] = connectEM.calculateSurfaceMyelinScore( axons, graph, borderMeta, heuristics ); % calculate myelin score for the dendrite class

disp('Done calculating TC/myelinated axons, creating isos')
%% export TC axons including myelin iso
Visualization.exportAggloToAmira(p, axons(isThalamocortical), fullfile(mainOutDir,'TC'), isoParamAx{:})
Visualization.exportAggloToAmira(p, myelinNeighSegments(isThalamocortical), fullfile(mainOutDir,'TCmyelin'), isoParamMy{:})


%% export myelinated axons including myelin iso, excluding tc axons
myAxNoTC = ~isThalamocortical & myelinFracAx > 0.1;
Visualization.exportAggloToAmira(p, axons(myAxNoTC), fullfile(mainOutDir,'AxMy_0.1b'), isoParamAx{:})
Visualization.exportAggloToAmira(p, myelinNeighSegments(myAxNoTC), fullfile(mainOutDir,'AxMymyelin_0.1b'), isoParamMy{:})
