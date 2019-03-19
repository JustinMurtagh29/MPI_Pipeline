% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
mix = struct('mean', {}, 'std', {}, 'coeff', {});
% Values for overall distribution of TC spine synapses. Used as bulk
mix(1).mean = -0.70; mix(1).std = 0.3; mix(1).coeff = 1;

pairCount = 5290;

info = Util.runInfo();
Util.showRunInfo(info);

%% Generate synapse pairs
clear cur*;
rng(0);

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

curName = tempname();
curInFile = strcat(curName, '.h5');
curOutFile = strcat(curName, '.out.h5');

for curIdx = 1:2
    curDset = sprintf('/log10Asi%d', curIdx);
    curData = pairT.areas(:, curIdx);
    
    h5create(curInFile, curDset, ...
        size(curData), 'Datatype', class(curData));
    h5write(curInFile, curDset, curData);
end

curCmd = fileparts(mfilename('fullpath'));
curCmd = fullfile(curCmd, 'clusterConnections.py');
curCmd = sprintf('%s "%s"', curCmd, curInFile);
[curExit, curResult] = system(curCmd, '-echo');

runT = table;
runT.means = transpose(h5read(curOutFile, '/mu'));
runT.stds = transpose(h5read(curOutFile, '/sigma'));
runT.coeffs = transpose(h5read(curOutFile, '/theta'));
