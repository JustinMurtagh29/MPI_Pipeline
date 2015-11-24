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
% agglo.borderMerged = agglomerateSG(graph, agglo.border, .95, 1, true);
% checkAggloConsitency(agglo);
% Lets agglomerate some nuclei (to somata) and seg. border IDs (to bigger neurites)
% agglo.borderGrown = agglomerateSG(graph, agglo.borderMerged, .95, 10, false);
% agglo.nucleiGrown = agglomerateSG(graph, agglo.nuclei, .95, 10, false);
load([p.saveFolder 'agglomeration' filesep 'nucleiVesselBorder.mat']);

% Seed 5 by 5 micron plane, see MH comment lablog
seedPlaneCenter = [4166 5376 2387];
expandVector = [floor(5000/11.24) 0 floor(5000/28)];
seedPlane(:,1) = seedPlaneCenter - expandVector;
seedPlane(:,2) = seedPlaneCenter + expandVector;
segIds = unique(readKnossosRoi(p.seg.root, p.seg.prefix, seedPlane, 'uint32', '', 'raw'));
segIds = segIds(segIds ~= 0);
segIds(ismember(segIds, cat(1,agglo.nuclei{:}, agglo.vessel{:}))) = [];
segIds = mat2cell(segIds, ones(length(segIds),1), 1);
% Agglomerate for 10 steps:
%aggloPlane = agglomerateSG(graph, segIds, .95, 10);
% Merge seeds which overlap in 1 segment
%uniqueAggloPlane = mergeSeeds(aggloPlane);
% Write results to skeleton structure and save
%skel = writeSkeleton(graph, aggloPlane, com);
%writeNml('/gaba/u/mberning/testPlane.nml', skel);
%skel = writeSkeleton(graph, uniqueAggloPlane, com);
%writeNml('/gaba/u/mberning/testPlaneMerged.nml', skel);

% Choose 100 random seeds from plane
idx = randi(length(segIds),100,1);
segIdsRand = segIds(idx);
p = [.9 .95 .98];
for p_i=1:length(p);
    aggloRand = agglomerateSG(graph, segIdsRand, p(p_i), 10);
    aggloRand = mergeSeeds(aggloRand);
    skel = writeSkeleton(graph, aggloRand, com);
    writeNml(['/gaba/u/mberning/testPlane' num2str(p(p_i), '%3.2f') '.nml'], skel);
end

aggloRand = agglomerateSG(graph, segIdsRand, .95, 100);
aggloRand = mergeSeeds(aggloRand);
skel = writeSkeleton(graph, aggloRand, com);
writeNml(['/gaba/u/mberning/testPlaneLong.nml'], skel);

aggloRand = agglomerateSG(graph, segIdsRand, .95, 10);
aggloRand = mergeSeeds(aggloRand);
skel = writeSkeleton(graph, aggloRand, com, true, segIds);
writeNml(['/gaba/u/mberning/testPlane' num2str(.95, '%3.2f') 'wQuerries.nml'], skel);
