function runChiasmaDetection(param)

dataDir = fullfile(param.saveFolder, 'aggloState/');

input = fullfile(dataDir, 'dendrites_flight_02.mat');

output = '/tmpscratch/scchr/L4/chiasmata/dendrites/';

% First round: 
% param.sphereRadiusOuter = 7000; % in nm
% param.sphereRadiusInner = 2000; % in nm
% param.minNodeDist = 4000; % in nm
% output = fullfile(output, 'round1/');
% Second round:
% param.sphereRadiusOuter = 20000; % in nm
% param.sphereRadiusInner = 2000; % in nm
% param.minNodeDist = 4000; % in nm
% output = fullfile(output, 'round2/');
% Third round:
% param.sphereRadiusOuter = 15000; % in nm
% param.sphereRadiusInner = 2000; % in nm
% param.minNodeDist = 4000; % in nm
% output = fullfile(output, 'round3/');
%Forth round:
param.sphereRadiusOuter = 10000; % in nm
param.sphereRadiusInner = 2000; % in nm
param.minNodeDist = 4000; % in nm
output = fullfile(output, 'round4/');

if ~exist(output, 'dir')
    mkdir(output)
end
    
job = connectEM.detectChiasmataSuperSuper( ...
    param, input, output);


%% Evaluation

aggloState = 'dendrites_flight_02';
% chiasmaId = '20171103T230345';
% chiasmaDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/round2';
% outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/round2';

% chiasmaId = '20171103T234130';
% chiasmaDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/round3';
% outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/round3';

chiasmaId = '20171104T105457';
chiasmaDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/round4';
outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/round4';

if ~exist(outputDir, 'dir')
    mkdir(outputDir)
end

tic;connectEM.evaluateChiasmataDetectionFunction(aggloState,chiasmaId,chiasmaDir,outputDir);toc;


%% Small subset using Marcels pseudorandom 50 agglos
outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/randomAgglos/round4';
aggloState = 'dendrites_flight_02';
Superagglos.get50RandDend(superagglos)
