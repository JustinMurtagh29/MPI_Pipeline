% filterPathLength 
% filters all agglomerates according to path length (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% agglos - agglomerates as sets of segments
% points - agglomerates as points from segment meta
% ulbound - 2x1 upperBound and lowerBound
% Outputs:
% filtered agglos & points, pathLengths (in um), rInd - indeces of filtered
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ filteredAgglos, filteredPoints, pathLengths, rInd ] = filterPathLength( agglos, points, ulbound, rInd )

% filter using path length
tic;
% create skeletons for contact filtered agglos
skeletons = [];
distthr = 15000;
for i=1:length(agglos)
   skel = skeleton();
   for j=1:length(agglos{i})
       skel.nodes = vertcat(skel.nodes, points{i}(j,:));
   end
   skel.nodes = vertcat(skel.nodes{:});
   skel.edges = Graph.getMST(bsxfun(@times, skel.nodes, [11.24 11.24 28]),distthr); 
   skeletons = vertcat(skeletons, skel);
end
disp('created skeletons for agglos.'); toc;

% calculate total path length per skeleton
pathLengths = zeros(length(skeletons),1);
for k=1:length(skeletons)
    pathLengths(k) = skeleton.physicalPathLength(skeletons(k).nodes, skeletons(k).edges, [11.24 11.24 28]);
end

% filter using upper lower bound
upperBound = ulbound(1);
lowerBound = ulbound(2);

% filtered set to be considered
filteredAgglos = vertcat(agglos(upperBound >= pathLengths & pathLengths >= lowerBound)); % with old file - two big apicals , agglos(62), agglos(75));
filteredPoints = vertcat(points(upperBound >= pathLengths & pathLengths >= lowerBound)); %, points(62), points(75));
rInd = rInd(upperBound >= pathLengths & pathLengths >= lowerBound);

f1 = length(filteredAgglos) / length(agglos);
fprintf('filtered agglos, fraction of total: %f (%d/%d)\n', f1, length(filteredAgglos), length(agglos));

disp('finished execution of filterPathLength.');
toc; disp('---------------------------------------');


end

