function dendritesNew = attachER(graph, dendrites, er)
% Attach ER that was FP of segment classifier back into dendrite class

erNeighbours = cellfun(@(x)unique(cat(2, graph.neighbours{x})), er, 'uni', 0);
% Find segment neighbouring each ER component for each dendrite component
for i=1:length(dendrites)
    dendriteNeighbours(i,:) = cellfun(@(x)sum(ismember(x, dendrites{i})), erNeighbours);
end
[maximumNeighbours, idx] = max(dendriteNeighbours, [], 1);
dendritesNew = dendrites;
for i=1:length(er)
    dendritesNew{idx(i)} = cat(1, dendritesNew{idx(i)}, er{i});
end

end

