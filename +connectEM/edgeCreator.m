function edges = edgeCreator(x,y)

    maxSize = max(numel(x),numel(y));
    maxSize = max(1, maxSize);
    edges = NaN(maxSize,2);
    edges(1:numel(x),1) = x;
    edges(1:numel(y),2) = y;

end

