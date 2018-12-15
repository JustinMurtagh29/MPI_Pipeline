function perm = linkageToPermutation(link)
    % perm = linkageToPermutation(link)
    %   linkageToPermutation(link) takes a hierarchical cluster tree (as
    %   produced by linkage) and returns a node permutation, such that
    %   similar nodes are close to each other. This function essentially
    %   performs the same operation as dendrogram, but without generating a
    %   plot and without being limited in the number of nodes.
    %
    % Written by
    %   Alessandro Motta <alessandro.motta@brain.mpg.de>
    maxId = link(end, 2) + 1;
    numNodes = maxId - size(link, 1);
    
    link(:, end) = numNodes + (1:size(link, 1));
    order = setdiff(transpose(link(:, 1:2)), link(:, end), 'stable');
    order = cat(1, flip(link(:, end)), order(:));
   [~, link] = ismember(link, order);
   
    S = repelem(link(:, end), 1, 2);
    T = link(:, 1:2);
    G = digraph(S, T);
    
    perm = dfsearch(G, 1);
    perm = order(perm);
    perm = perm(perm <= numNodes);
end
