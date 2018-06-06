% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

minSynPost = 10;
info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);
conn = connectEM.Connectome.prepareForSpecificityAnalysis(conn);
axonClasses = unique(conn.axonMeta.axonClass);

%% build class connectome
classConnectome = ...
    connectEM.Connectome.buildClassConnectome( ...
        conn, 'targetClasses', [], 'axonClasses', axonClasses);

[dendMeta, classConnectome] = ...
    connectEM.Connectome.prepareForFullCellInputAnalysis( ...
        conn.denMeta, classConnectome);

classConnectome = transpose(classConnectome);
axonClasses = reshape(axonClasses, 1, []);

%% build dendrite class(es)
dendClasses = struct;
dendClasses(1).ids = find( ...
    dendMeta.targetClass ~= 'Somata' ...
  & dendMeta.targetClass ~= 'AxonInitialSegment' ...
  & dendMeta.targetClass ~= 'FullInput' ...
  & dendMeta.synCount >= minSynPost);
dendClasses(1).nullIds = dendClasses(1).ids;
dendClasses(1).title = sprintf( ...
    'Dendrites with ≥ %d synapses (n = %d)', ...
    minSynPost, numel(dendClasses(1).ids));

dendClasses(2).ids = find( ...
    dendMeta.targetClass == 'FullInput' ...
  & dendMeta.synCount >= 500);
dendClasses(2).nullIds = dendClasses(2).ids;
dendClasses(2).title = sprintf( ...
    'Whole cells with ≥ %d synapses (n = %d)', ...
    500, numel(dendClasses(2).ids));
    
%% plot
for curIdx = 1:numel(dendClasses)
    connectEM.Specificity.plotRatioDist( ...
        info, classConnectome, axonClasses, dendClasses(curIdx));
end
