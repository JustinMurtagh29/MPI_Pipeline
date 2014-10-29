function g = constructGraph_v2(raw, aff, seg)

% Initalize structure
g = mygraph();

% Find nearest neighbours in segmentation and determine weights according to classification
g = calcNeighbourList_v2(g, raw, aff, seg);

end
