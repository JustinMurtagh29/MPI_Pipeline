function t = thresholdsForFocusedAnnotationScoreBased( ...
    scores, target, reqP, reqR)
%THRESHOLDSFORFOCUSEDANNOTATIONSCOREBASED Calculate the thresholds for
%focused annotation using the test set prediction scores.
% INPUT scores: [Nx2] or [Nx1] float
%           The test set scores
%       target: [Nx1] int or logical
%           The test set labels. In case those are int then the targets are
%           interpreted as groups and only the maximal score for each group
%           is kept.
%       reqP: double
%           Required precision. If > 1 then it is normalized by 100.
%       reqR: double
%           Required recall. If > 1 then it is normalized by 100.
% OUTPUT t: [1x2] double
%           Lower and the upper threshold for focused annotation.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if size(scores, 2) == 2
    scores = max(scores, [], 2);
end

%normalize RP
if reqP > 1; reqP = reqP /100; end
if reqR > 1; reqR = reqR /100; end

if ~islogical(target)
    group = Util.group2Cell(target);
    
    %delete zero group
    group(1) = [];
    
    %delete instance with maximal score from group
    for i = 1:length(group)
        [~,idx] = max(scores(group{i}));
        group{i}(idx) = [];
    end
    
    %delete all remaining elements in groups from scores and targets
    toDel = cell2mat(group(~cellfun(@isempty,group)));
    scores(toDel) = [];
    target(toDel) = [];
    target = target > 0;
end

%determine R threshold (lower)
[scores, sI] = sort(scores, 'ascend');
target = target(sI);
r = cumsum(target, 'reverse')./sum(target);
lowerTIdx = find(r >= reqR, 1, 'last');

%determint P threshold (upper)
%this assumes that all interfaces between the lower and the upper threshold
%are correct and calculates precision including those correct ones
for i = lowerTIdx:length(scores)
    p = sum(target(lowerTIdx:end))/(sum(1 - target(i:end)) + sum(target(lowerTIdx:end)));
    if p > reqP
        break
    end
end
upperTIdx = i;

% %determine P threshold (upper)
% p = cumsum(target, 'reverse')./(length(target):-1:1)';
% upperTIdx = find(p > reqP, 1, 'first');

t = scores([lowerTIdx, upperTIdx]);
end

