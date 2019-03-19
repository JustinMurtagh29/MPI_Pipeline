% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
mix = struct('mean', {}, 'std', {}, 'coeff', {});
% Values for overall distribution of TC spine synapses. Used as bulk
mix(1).mean = -0.70; mix(1).std = 0.3; mix(1).coeff = 0.8;
mix(2).mean = -0.25; mix(2).std = 0.2; mix(2).coeff = 0.2;

pairCount = 5290;

info = Util.runInfo();
Util.showRunInfo(info);

%% Generate synapse pairs
clear cur*;

curClassEdges = [mix.coeff];
curClassEdges = curClassEdges / sum(curClassEdges);
curClassEdges = cumsum([0; curClassEdges(:)]);

pairT = table;
pairT.classId = reshape(linspace(0, 1, pairCount), [], 1);
pairT.classId = discretize(pairT.classId, curClassEdges);

pairT.areas = randn(height(pairT), 2);
pairT.areas = cat(1, mix(pairT.classId).std) .* pairT.areas;
pairT.areas = cat(1, mix(pairT.classId).mean) + pairT.areas;

%% Run inference
clear cur*;

curH5File = strcat(tempname, '.h5');
for curIdx = 1:2
    curDset = sprintf('/log10Asi%d', curIdx);
    curData = pairT.areas(:, curIdx);
    
    h5create(curH5File, curDset, size(curData), 'Datatype', class(curData));
    h5write(curH5File, curDset, curData);
end

curCmd = fileparts(mfilename('fullpath'));
curCmd = fullfile(curCmd, 'clusterConnections.py');
curCmd = sprintf('%s "%s"', curCmd, curH5File);
[curExit, curResult] = system(curCmd, '-echo');
