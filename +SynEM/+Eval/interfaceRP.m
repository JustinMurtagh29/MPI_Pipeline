function [rp, thresholds, auc] = interfaceRP( target, scores, thresholds )
%INTERFACERP RP curve of two-class classification problem with optional
%grouping.
% INPUT target: [Nx1] Logical or int
%           Target labels (see also Util.confusionMatrix) or grouped
%           targets.
%       scores: [Nx1] float of prediction scores.
%           Prediction scores.
%       thresholds: (Optional) [Nx1] double
%           The dection thresholds for which rp values are evaluated. If
%           this input is provided the rows in rp are not uniquified.
%           (Default: min(scores):0.01:max(scores))
% OUTPUT rp: [Nx2] double
%           Pairs of recall precision values in each row. Rows are
%           uniquified if no threshold is given as input.
%        thresholds: [Nx1] double 
%           The thresholds for the corresponding rows in rp.
%       auc: double
%           Area under curve determined using trapezoidal numerical
%           integration.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if isrow(scores)
    scores = scores';
end

if ~exist('thresholds', 'var') || isempty(thresholds)
    thresholds = unique(scores);
    thresholds = (min(thresholds):0.01:max(thresholds))';
    uniquify = true;
else
    uniquify = false;
end

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

if uniquify
    [rp, idx] = unique(rp, 'rows', 'stable');
    thresholds = thresholds(idx);
end
auc = trapz(rp(end:-1:1,1), rp(end:-1:1,2));
end

