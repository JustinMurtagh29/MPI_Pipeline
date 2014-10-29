function g = constructGraph(raw, aff, seg)

% Initalize structure
g = mygraph();

% Find nearest neighbours in segmentation and determine weights according to classification
g = calcNeighbourList(g, raw, aff, seg);

end
