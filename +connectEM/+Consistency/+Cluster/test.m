% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
mix = struct('mean', {}, 'std', {}, 'coeff', {});

i = 1;
% mix(i).mean = -1.00; mix(i).std = 0.2; mix(i).coeff = 0.1; i = i + 1; %#ok
mix(i).mean = -0.70; mix(i).std = 0.3; mix(i).coeff = 0.8; i = i + 1; %#ok
mix(i).mean = -0.25; mix(i).std = 0.2; mix(i).coeff = 0.2; i = i + 1; %#ok
clear i;

pairCount = 5290;

info = Util.runInfo();
Util.showRunInfo(info);

%% Generate synapse pairs
clear cur*;
rng(0);

pairT = table;
[pairT.areas, pairT.classId] = ...
    connectEM.Consistency.Simulation.sampleGmm( ...
        mix, 'numSynapses', 2, 'numConnections', pairCount);

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
