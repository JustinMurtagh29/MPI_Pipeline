function agglomerationPostHocTwo(todo, graph, gridAgglo_05, segmentMeta, borderMeta, todo2)

    if ~exist('todo2', 'var')
        todo2 = find(ismember(graph.edges, todo, 'rows'));
    end
    todo2(isnan(graph.borderIdx(todo2))) = [];
    graph.prob(todo2) = 42; %magic number if we need to find this later
    segmentMeta.voxelCount(flatten(graph.edges(todo2,:))) = Inf;
    borderMeta.borderSize(graph.borderIdx(todo2)) = Inf;
    connectEM.agglomerationModify(gridAgglo_05{564}, ...
        '/gaba/scratch/kboerg/aggloSearch/aggloPostHoc.mat', ...
        graph, ...
        struct('segmentMeta', segmentMeta, 'borderMeta', borderMeta);
end
function y =flatten(x)
    y=x(:);
end
