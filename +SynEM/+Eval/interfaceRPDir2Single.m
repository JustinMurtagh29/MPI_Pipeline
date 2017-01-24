function rp = interfaceRPDir2Single( target, scores, mode )
%INTERFACERPDIR2SINGLE RP values for directed prediction but evaluated only
%whether an interface was corrently classified independent of its
%direction.
% INPUT target: [Nx1] cell array
%           Each cell contains the labels for one cube.
%       scores: [Nx1] double
%           Prediction scores for all cubes.
%       mode: (Optional) string
%           Input mode of target. Options are
%           'single': (Default) Target labels contains one label for each
%               interface.
%           'directed': Target labels are for direction labels
% OUTPUT rp: [Nx2] double
%           Pairs of recall precision values in each row.
%        thresholds: [Nx1] double 
%           The thresholds for the corresponding rows in rp.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('mode','var') || isempty(mode)
    mode = 'single';
end

l = cellfun(@length,target);
switch mode
    case 'single'
        scores = mat2cell(scores, 2.*l, 1);
        scores = cellfun(@(x)reshape(x,[],2), scores, ...
            'UniformOutput', false);
    case 'directed'
        scores = mat2cell(scores, l, 1);
        scores = cellfun(@(x)reshape(x,[],2), scores, ...
            'UniformOutput', false);
        target = cellfun(@(x)any(reshape(x,[],2),2), target, ...
            'UniformOutput', false);
    otherwise
        error('Unknown mode %s.', mode);
end
target = cell2mat(target);
scoresM = cell2mat(scores);

thresholds = unique(scoresM(:));
thresholds = (min(thresholds):0.01:max(thresholds))';
targetScores = max(scoresM(target,:),[],2);
%add target scores to thresholds
thresholds = unique([thresholds; targetScores]); 
rp = zeros(length(thresholds),2);
for i = 1:length(thresholds)
    pred = cell2mat(cellfun(@(x)any(x >= thresholds(i),2), scores, ...
        'UniformOutput', false));
    [~,p, r] = SynEM.Util.confusionMatrix(target, pred);
    rp(i,:) = [r, p];
end

end

