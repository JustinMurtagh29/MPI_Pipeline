function [rp, thresholds] = interfaceRP( target, scores )
%INTERFACERP RP curve of two-class classification problem with optional
%grouping.
% INPUT target: [Nx1] Logical or int
%           Target labels (see also Util.confusionMatrix) or grouped
%           targets.
%       scores: [Nx1] float of prediction scores.
%           Prediction scores.
% OUTPUT rp: [Nx2] double
%           Pairs of recall precision values in each row.
%        thresholds: [Nx1] double 
%           The thresholds for the corresponding rows in rp.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isrow(scores)
    scores = scores';
end

thresholds = unique(scores);
thresholds = (min(thresholds):0.01:max(thresholds))';

%add target scores to thresholds
if islogical(target)
    targetScores = scores(target);
    thresholds = unique([thresholds; targetScores]);
end

rp = zeros(length(thresholds),2);
for i = 1:length(thresholds)
    y = scores >= thresholds(i);
    [~,p,r] = SynEM.Util.confusionMatrix(target, y );
    rp(i,:) = [r, p];
end


end

