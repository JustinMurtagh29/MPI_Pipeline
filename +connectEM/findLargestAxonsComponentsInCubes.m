
% Needs graph and segmentMeta with transposed .point
% Note: 00023 is 90% 00024 is 93%
a = load('/gaba/scratch/mberning/aggloSearch/00024.mat');                                                                                
axonsSize = cellfun(@(x)sum(segmentMeta.voxelCount(x)), a.axonsFinal);                                                                      
[axonSize, idx] = sort(axonsSize, 'descend');
axonsFinal = a.axonsFinal(idx);
Superagglos.skeletonFromAgglo(graph.edges, segmentMeta, axonsFinal(1:200), 'axonsLargest93_', '/gaba/scratch/mberning/alignmentSkeletons/'); 
save('/gaba/u/mberning/repos/pipeline/+connectEM/evaluationData/gliaER/cmAndErAnnotations.mat', 'axonsFinal');

