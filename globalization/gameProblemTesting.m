function gameProblemTesting(p, coord, lowerT, upperT)
% Some testing for annotating single axon or dendrite
% Input: Usual parameter structure, coordinate for start, lower and upper threshold

% Look up segment ID at inital location
segId = readKnossosRoi(p.seg.root, p.seg.prefix, [coord; coord])

% Load supervoxel in global IDs (edges & prob before joining, (edge/prob)Remaining after GP and correspondence application)
load([p.saveFolder 'graphNew.mat']);

% Write skeletons (later here, game problems) according to pruned graph


end
