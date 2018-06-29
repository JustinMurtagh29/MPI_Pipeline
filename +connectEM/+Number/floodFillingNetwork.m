% All numbers based on
%   https://arxiv.org/abs/1611.00421v1
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

info = Util.runInfo();
Util.showRunInfo(info);

%% Configuration
pathLenMm = 97;

axonFrac = 0.34;
axonRecall = 0.68;

dendFrac = 1 - axonFrac;
dendRecall = 0.48;

mergeErrors = 3E0;
splitErrors = 3E3;

%% Error rate
axonPathLenMm = axonFrac * pathLenMm %#ok
dendPathLenMm = dendFrac * pathLenMm %#ok

reconAxonPathLenMm = axonRecall * axonPathLenMm %#ok
reconDendPathLenMm = dendRecall * dendPathLenMm %#ok

reconPathLen = reconAxonPathLenMm + reconDendPathLenMm %#ok
errorsPerMm = (mergeErrors + splitErrors) / reconPathLen %#ok
