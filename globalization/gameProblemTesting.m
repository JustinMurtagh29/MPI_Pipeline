function gameProblemTesting(p, coord, lowerT, upperT)
% Some testing for annotating single axon or dendrite
% Input: Usual parameter structure, coordinate for start, lower and upper threshold

% Look up segment ID at inital location
segId = readKnossosRoi(p.seg.root, p.seg.prefix, [coord; coord])

% Load correspondences in global IDs, e.g. how to continue over cube borders (variable components)


% Write skeletons (later here, game problems) according to pruned graph


end
