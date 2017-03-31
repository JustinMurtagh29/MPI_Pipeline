function [graphCut, mapping] = cutGraph(p, segmentMeta, borderMeta, outputFolder, borderSizeThreshold, segmentSizeThreshold)
% Restrict graph based on heuristics results

% Segments classified by heuristics
load([p.saveFolder 'heuristicResult.mat']);
assert(length(segIds) == max(segIds));
vesselIds = vesselScore > 0.5;
nucleiIds = nucleiScore > 0.5 & ~vesselIds;
myelinIds = myelinScore > 0.5 & ~vesselIds & ~nucleiIds;
heuristicIds = vesselIds | nucleiIds | myelinIds;
% Keep only segments larger than segmentSizeThreshold voxel and ...
% that have more than 50% connection probability to another segment after removing border smaller than borderSizeThreshold
load([p.saveFolder  'globalGPProbList.mat']);
load([p.saveFolder 'globalEdges.mat']);
edgeIdx = find(borderMeta.borderSize > borderSizeThreshold);
remainingEdges = edges(edgeIdx, :);
remainingProb = prob(edgeIdx);
% Add correspondences
corrEdges = Seg.Global.getGlobalCorrespondences(p);
corrProb  = ones(size(corrEdges, 1), 1);
remainingEdges = [remainingEdges; corrEdges];
remainingProb  = [remainingProb ; corrProb];
maxProb = accumarray(cat(1,remainingEdges(:,1),remainingEdges(:,2)), cat(1,remainingProb, remainingProb),[segmentMeta.maxSegId 1], @max);
smallIds = segmentMeta.voxelCount <= segmentSizeThreshold & ~heuristicIds;
lowProbIds = segmentMeta.voxelCount > segmentSizeThreshold & maxProb <= 0.5 & ~heuristicIds;
% Write mapping script for visualization
mapping = {vesselIds nucleiIds myelinIds smallIds lowProbIds};
mapping = cellfun(@find, mapping, 'uni', 0);
mappingFile = [outputFolder 'cutFromGraph_' num2str(borderSizeThreshold, '%.5i') '_' num2str(segmentSizeThreshold, '%.5i') '.txt'];
script = WK.makeMappingScript(segmentMeta.maxSegId, mapping, false);
fileHandle = fopen(mappingFile, 'w');
fwrite(fileHandle, script);
fclose(fileHandle);
% Remove heuristics and small or disconnected segments from graph for now
removedIds = cat(1, mapping{:});
keptIds = setdiff(1:segmentMeta.maxSegId, removedIds);
keepEdgeIdx = all(ismember(remainingEdges, keptIds), 2);
graphCut.edges = remainingEdges(keepEdgeIdx,:);
graphCut.prob = remainingProb(keepEdgeIdx);
% Sort edges again (e.g. make correspondences find their place
[graphCut.edges, edgeRows] = sortrows(graphCut.edges);
graphCut.prob = graphCut.prob(edgeRows);

end

