function [ fps, fns ] = getPotentialClassificationErrors( y, scores, mode )
%GETPOTENTIALCLASSIFICATIONERRORS Get class instances with highest score in
% in negative class (potential fps) and with lowest score in positive class
% (potential fns).
% INPUT y: [Nx1] logical array of ground truth labels.
%       scores: [Nx1] numerical array of prediction scores.
%       mode: (Optional) See getTrainingDataFrom.
%           (Default: 'direction')
% OUTPUT fps: Table containing the index of the potential fp in y and their
%           scores.
%        fns: Table containing the index of the potential fn in y and their
%           scores.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('mode','var') || isempty(mode)
    mode = 'direction';
end
%reshape to p
switch mode
    case 'direction'
        y = any(reshape(y,length(y)/2,2),2);
        scores = max(reshape(scores,length(scores)/2,2),[],2);
    otherwise
        error('Mode %s not implemented', mode);
end

%fps
[sScores, I] = sort(scores,'descend');
fps = [I(~y(I)),sScores(~y(I))];

%fns
[sScores, I] = sort(scores,'ascend');
fns = [I(y(I)),sScores(y(I))];

%convert to table
fps = array2table(fps,'VariableNames',{'Index','Score'});
fns = array2table(fns,'VariableNames',{'Index','Score'});


end

