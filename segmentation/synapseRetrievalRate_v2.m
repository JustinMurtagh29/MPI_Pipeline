function synRate = synapseRetrievalRate_v2( seg, skeletons )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


totalSkeletons = 0;
skeletonNodes = zeros(length(skeletons),1);
segObjectAtNode = cell(length(skeletons),1);
for s=1:length(skeletons)
    segObjectAtNode{s} = seg(skeletons{s}.nodes(:,1:3));
    skeletonNodes(s) = size(skeletons{s}.nodes,1);
    totalSkeletons = totalSkeletons + 1;
end
totalNodes = sum(skeletonNodes);

truePos2 = 0;
truePosCounter = 0;
falsePos2 = 0;
falsePosCounter = 0;
% Brute Force Methode (iterate over all node pairs)
for s1=1:length(skeletons)
    for s2=1:length(skeletons)
        for obj1=1:length(segObjectAtNode{s1})
            for obj2=1:length(segObjectAtNode{s2})
                if segObjectAtNode{s1}(obj1) ~= 0
                    if s1 == s2 && obj1 ~= obj2
                        truePosCounter = truePosCounter + 1;
                        if segObjectAtNode{s1}(obj1) == segObjectAtNode{s2}(obj2)
                            truePos2 = truePos2 + 1;
                        end
                    end
                    if s1 ~= s2
                        falsePosCounter = falsePosCounter + 1;
                        if segObjectAtNode{s1}(obj1) == segObjectAtNode{s2}(obj2)
                            falsePos2 = falsePos2 + 1;
                        end
                    end
                end
            end
        end
    end
end
truePos1 = truePos2 ./ truePosCounter;
falsePos1 = falsePos2 ./ falsePosCounter;

% TO DO: Correct for BBox and Neuron Size effects

% Calculate synapse retreival rate
npRate = (fp1all + totalNodes/totalSkeletons - truePosAllCorrected) / (totalNodes/totalSkeletons);
npAcc = 1 - npRate;
cAcc = npAcc ^ 2;
synRate = 1 - (1 - cAcc) / totalNodes;

end

