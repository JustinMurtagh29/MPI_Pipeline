function [edges, idxMultipleHits] = edgeCreator(x,y)

    idxMultipleHits = numel(y) > 1;
    if idxMultipleHits
        % Ignore multiple ending queries for now
        edges = zeros(0,2);
    else
        edges = NaN(1,2);
        edges(1:numel(x),1) = x;
        edges(1:numel(y),2) = y;
    end

end

