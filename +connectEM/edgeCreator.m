function [edges, idxMultipleHits] = edgeCreator(sources, targets, allowMultiEnd)
    % [edges, idxMultipleHits] = edgeCreator(sources, targets)
    %   Returns the directed edges induced by the given `sources` and
    %   `targets`.
    %
    %   If `targets` contains multiple values it is considered to be
    %   invalid. In this case the `idxMultipleHits` flag is `true` and the
    %   `edges` list is empty (i.e., has zero rows). This is not the case
    %   if the allowMultiEnd is set to true
    %
    %   In practice `edges` contains at most one edge because `sources`
    %   also should contain at most one value.
    %
    % Written by
    %   Christian Schramm <christian.schramm@brain.mpg.de>
    %   Manuel Berning <manuel.berning@brain.mpg.de>
    
    if ~exist('allowMultiEnd','var') || isempty(allowMultiEnd)
        allowMultiEnd = 0;
    end
    % NOTE(amotta): Feel free to remove this line, Christian and Manuel.
    assert(numel(sources) < 2);
    
    % Check whether multiple endings were reached.
    idxMultipleHits = numel(targets) > 1;
 
    if idxMultipleHits
        if allowMultiEnd
            edges = [repmat(sources,numel(targets),1) , targets(:)];
        else
            % Flight paths which reach multiple ends (i.e., agglomerates or
            % endings) are ignored. Return empty edge list instead.
            edges = zeros(0,2);
        end
    else
        edges = NaN(1,2);
        edges(1:numel(sources),1) = sources;
        edges(1:numel(targets),2) = targets;
    end

end

