function [ rp, t ] = binaryRP( target, scores )
%BINARYRP Precision and recall for a binary classification task.
% INPUT target: [Nx1] logical or int
%           Binary target labels. If target is an integer than it is
%           assumed that each integer defines a grouping and only one
%           member out of a group needs to be found in order to consider
%           the whole group detected (corresponding to consider only the
%           largest score of a group and deleteting all others).
%       scores: [Nx1] float of prediction scores.
%           Prediction scores.
% OUTPUT rp: [Nx2] double
%           Pairs of recall and precision values in each row.
%        t: [Nx1] double
%           Score threshold for the corresponding rp pair.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

% keep only largest element form each group
if ~islogical(target)
    group = Util.group2Cell(target);
    
    % delete zero group
    group(1) = [];
    
    % delete instance with maximal score from group
    for i = 1:length(group)
        [~,idx] = max(scores(group{i}));
        group{i}(idx) = [];
    end
    
    % delete all remaining elements in groups from scores and targets
    toDel = cell2mat(group(~cellfun(@isempty,group)));
    scores(toDel) = [];
    target(toDel) = [];
    target = target > 0;
end

rp = zeros(length(scores), 2);

% recall
[t, sI] = sort(scores, 'ascend');
target = target(sI);
rp(:,1) = cumsum(target, 'reverse')./sum(target);

% precision
rp(:,2) = cumsum(target, 'reverse')./(length(target):-1:1)';
end

