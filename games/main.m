% Script to keep in mind what I have done for game problem generation

% Load parameter file from new restricted version of pipeline
load('/gaba/u/mberning/results/pipeline/20151111T183414/allParameter.mat');

% Calculate global graph
% collectGlobalGraph(p);
% Just load result now
load([p.saveFolder 'graph.mat']);

% Calculate CoM of all segments in supervoxel graph
% maxID = max(graph.edges(:));
% collectGlobalCoMs(p, maxID);
% Just load result now
load([p.saveFolder 'CoM.mat']);

% Commented out as stored result is used
% agglo.nuclei = determineSegIdsFromMag4(p, p.nuclei); 
% agglo.vessel = determineSegIdsFromMag4(p, p.vessel);
% excludedID = cat(1,agglo.nuclei{:}, agglo.vessel{:});
% agglo.border = determineSegIdsAtSGborders(p, excludedID);
agglo.borderMerged = agglomerateSG(graph, agglo.border, .95, 1, true);
load([p.saveFolder 'agglomeration' filesep 'nucleiVesselBorder.mat']);
agglo = rmfield(agglo, 'border');
checkAggloConsitency(agglo);

% Lets agglomerate some nuclei (to somata) and seg. border IDs (to bigger neurites)
agglo.borderGrown = agglomerateSG(graph, agglo.borderMerged, .95, 10, false);
agglo.nucleiGrown = agglomerateSG(graph, agglo.nuclei, .95, 10, false);


