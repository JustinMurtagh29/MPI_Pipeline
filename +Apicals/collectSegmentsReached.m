% collectSegmentsReached 
% collect the segments in a specific dendrite agglo have been (if at all)
% reached by a certain spine head (designed for L4)
% -------------------------------------------------------------------------
% Inputs:
% spineEdges - edges from spine heads attachment by Alessandro
% shAgglos - spine head agglos from spine heads attachment by Alessandro
% Outputs:
% segmentsReached - spine head to dend agglo segment (if attached)
% -------------------------------------------------------------------------
% Author: Matej Zecevic <matej.zecevic@brain.mpg.de>
function [ segmentsReached ] = collectSegmentsReached( spineEdges, shAgglos )
tic;
% use Alessandros edge list to get the segments each spinehead encounters
segmentsReached = zeros(length(spineEdges),1);
for i=1:length(segmentsReached)
    ind = find(spineEdges{i} == 0);
    if length(ind) == 20 % if all are 0
       segmentsReached(i) = shAgglos{i}(1);
    end
    if isempty(ind) % if non is 0
        segmentsReached(i) = spineEdges{i}(end);
    else
        for j=1:10 % find last element
            if spineEdges{i}(j,2) == 0 && spineEdges{i}(j,1) ~= 0
                segmentsReached(i) = spineEdges{i}(j,1);
            elseif spineEdges{i}(j,2) ~= 0 && spineEdges{i}(j+1,1) == 0
                segmentsReached(i) = spineEdges{i}(j,2);
            end
        end
    end
%     if mod(i, 1000) == 0
%         fprintf('finished %d \n', i);
%     end
end

disp('finished execution of collectSegmentsReached.');
toc; disp('---------------------------------------');

end

