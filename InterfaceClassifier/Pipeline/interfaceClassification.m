function [X,  interfaces ] = interfaceClassification( raw, seg, edges, borders, featureMap, classifier )
%INTERFACECLASSIFICATION Calculation of interfaces, features and subsequent
% classifiation.
% INPUT raw: 3d array of single containing the raw data.
%       seg: 3d array of integer containing a segmentation.
%       edges: [Nx2] array of integer containing the edges in the adjacency
%           graph contained in seg.
%       borders: Struct containing the borders between edges.
%       featureMap: see makeFeatureMap
%       classifier: Interface classifier.
% OUTPUT scores: [Nx2] array where N = size(edges,1) containing the synapse
%           score for each edge in edges. If the edge size is below the
%           area threshold than NaN is saved.
%        X: [NxM] feature matrix where N = 2*size(edges,1). Rows correspond
%           to interfaces, columns to features.
%        interfaces: see calculateInterfaces.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

[interfaces, intIdx] = calculateInterfaces(seg, edges, borders, featureMap.areaThreshold, featureMap.voxelSize, featureMap.rinclude);
X = calculateFeaturesOld(raw, interfaces, featureMap);

end
