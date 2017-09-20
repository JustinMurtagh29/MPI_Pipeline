function [edges, idxMultipleHits] = edgeCreator(sources, targets)
    % [edges, idxMultipleHits] = edgeCreator(sources, targets)
    %   Returns the directed edges induced by the given `sources` and
    %   `targets`.
    %
    %   If `targets` contains multiple values it is considered to be
    %   invalid. In this case the `idxMultipleHits` flag is `true` and the
    %   `edges` list is empty (i.e., has zero rows).
    %
    %   In practice `edges` contains at most one edge because `sources`
    %   also should contain at most one value.
    %
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    
    % NOTE(amotta): Feel free to remove this line, Christian and Manuel.
    assert(numel(sources) < 2);
    
    % Check whether multiple endings were reached.
    idxMultipleHits = numel(targets) > 1;
    
    if idxMultipleHits
        % Flight paths which reach multiple ends (i.e., agglomerates or
        % endings) are ignored. Return empty edge list instead.
        edges = zeros(0,2);
    else
        edges = NaN(1,2);
        edges(1:numel(sources),1) = sources;
        edges(1:numel(targets),2) = targets;
    end

end

