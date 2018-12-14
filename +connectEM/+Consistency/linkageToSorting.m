function sorting = linkageToSorting(links)
    maxId = links(end, 2) + 1;
    numNodes = maxId - size(links, 1);
    
    links(:, end) = numNodes + (1:size(links, 1));
    order = setdiff(links(:, 1:2)', (numNodes + 1):maxId, 'stable');
    order = cat(1, rot90((numNodes + 1):maxId), order(:));
   [~, links] = ismember(links, order);
   
    S = repelem(links(:, end), 1, 2);
    T = links(:, 1:2);
    G = digraph(S, T);
    
    sorting = dfsearch(G, 1);
    sorting = order(sorting);
    sorting = sorting(sorting <= numNodes);
end
