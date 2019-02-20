% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
%   Sahil Loomba <sahil.loomba@brain.mpg.de>
% do Hierarchical clustering per dataset cube and write graph file with edges and scores
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

margin = param.tileBorder(:, 2);
assert(isequal(margin, [256; 256; 128]));

info = Util.runInfo();
Util.showRunInfo(info);

%% Run per-cube agglomeration
taskArgs = arrayfun( ...
    @(local) {local.bboxSmall, fullfile(local.saveFolder, 'graphH.mat')}, ...
    param.local, 'UniformOutput', false);
sharedArgs = {param, margin};
sharedArgLocs = [1, 4];

job = Cluster.startJob( ...
    @jobWrapper, ...
    taskArgs, 'numOutputs', 0, ...
    'sharedInputs', sharedArgs, ...
    'sharedInputsLocation', sharedArgLocs, ...
    'cluster', {'priority', 100, 'time', '24:00:00', 'memory', 48});
Cluster.waitForJob(job);

function jobWrapper(param, bbox, graphFile, margin)

    [edges, scores] = connectEM.Hierarchical.runBox(param, bbox, 'margin', margin);

    %% Build graph
    [curCount, ~] = cellfun(@size, edges);
    scores = repelem(scores, curCount, 1);
    edges = cat(1, edges{:});
    
    graph = table;
    graph.edge = edges;
    graph.score = scores;
    clear edges scores;
    
    graph = sortrows(graph, 'score', 'descend');
    [~, curUni] = unique(graph.edge, 'rows', 'stable');
    graph = graph(curUni, :);
    
    Util.save(graphFile, graph);
end
