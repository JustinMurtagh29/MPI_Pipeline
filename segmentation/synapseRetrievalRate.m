function synRate = synapseRetrievalRate( seg, skeletons )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


totalSkeletons = 0;
skeletonNodes = zeros(length(skeletons),1);
segObjectAtNode = cell(length(skeletons),1);
for s=1:length(skeletons)
    if ~isempty(skeletons{s}.nodes)
        segObjectAtNode{s} = seg(skeletons{s}.nodes(:,1:3));
        skeletonNodes(s) = size(skeletons{s}.nodes,1);
        totalSkeletons = totalSkeletons + 1;
    end
end
totalNodes = sum(skeletonNodes);

truePos1 = zeros(length(skeletons),1);
truePos2 = zeros(length(skeletons),1);
falsePos2 = zeros(length(skeletons),1);
for s1=1:length(segObjectAtNode)
    for obj=1:length(segObjectAtNode{s1})
        if segObjectAtNode{s1}(obj) ~= 0
            % Count all objects that are the same as segObjectAtNode obj &
            % correct for not counting the same object at same node as TP
            truePos2(s1) = truePos2(s1) + sum(segObjectAtNode{s1}(obj) == segObjectAtNode{s1}) - 1;
        end
    end
    % Normalize to total number of nodes in skeleton
    truePos1(s1) = truePos2(s1) ./ skeletonNodes(s1);
    for obj=1:length(segObjectAtNode{s1})
        for s2=1:length(segObjectAtNode)
            if segObjectAtNode{s1}(obj) ~= 0 && s1 ~= s2
                % Count all objects that are the same as segObjectAtNode obj
                falsePos2(s1) = falsePos2(s1) + sum(segObjectAtNode{s1}(obj) == segObjectAtNode{s2});
            end
        end
    end
    
end
% Mean (of truePos of a skeleton weighted by node number)
truePosAll = sum(truePos1 .* skeletonNodes) ./ totalNodes;
% Correct for bounding box & neuron size effects
truePosAllCorrected = 1 ./ ((1 ./ truePosAll) - (totalSkeletons ./ totalNodes));
% Normalize to total number of ALL nodes
falsePosAll = sum(falsePos2) ./ totalNodes;

% Calculate synapse retreival rate
npRate = (falsePosAll + totalNodes/totalSkeletons - truePosAllCorrected) / (totalNodes/totalSkeletons);
npAcc = 1 - npRate;
cAcc = npAcc ^ 2;
synRate = 1 - (1 - cAcc) / totalNodes;

end

