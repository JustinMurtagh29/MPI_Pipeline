function runChiasmaDetection(param)

dataDir = fullfile(param.saveFolder, 'aggloState/');

input = fullfile(dataDir, 'dendrites_flight_02.mat');

output = '/tmpscratch/scchr/L4/chiasmata/dendrites/';

% First round: p.sphereRadiusOuter = 7000; % in nm
%              p.sphereRadiusInner = 2000; % in nm
%              p.minNodeDist = 4000; % in nm
% output = fullfile(output, 'round1/');
% Second round:
% param.sphereRadiusOuter = 20000; % in nm
% param.sphereRadiusInner = 2000; % in nm
% param.minNodeDist = 4000; % in nm
% output = fullfile(output, 'round2/');
% Third round:
p.sphereRadiusOuter = 15000; % in nm
p.sphereRadiusInner = 2000; % in nm
p.minNodeDist = 4000; % in nm
output = fullfile(output, 'round3/');



if ~exist(output, 'dir')
    mkdir(output)
end
    
job = connectEM.detectChiasmataSuperSuper( ...
    param, input, output);


%% Evaluation

aggloState = fullfile(dataDir,'dendrites_flight_02.mat');
chiasmaId = '20171103T230345';
outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/round2';
% outputDir = '/tmpscratch/scchr/L4/chiasmata/dendrites/evaluation/round3';
% chiasmaId = '20171103T230625';

tic;connectEM.evaluateChiasmataDetectionFunction(aggloState,chiasmId,chiasmDir,outputDir);toc;
