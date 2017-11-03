function runChiasmaDetection(param)

dataDir = fullfile(param.saveFolder, 'aggloState/');

input = fullfile(dataDir, 'dendrites_flight_02.mat');

output = fullfile(dataDir, 'chiasmata/dendrites');

% First round: p.sphereRadiusOuter = 7000; % in nm
%              p.sphereRadiusInner = 2000; % in nm
output = fullfile(output, 'round1/');

if ~exist(output, 'dir')
    mkdir(output)
end
    
job = detectChiasmataSuperSuper( ...
    p, input, output);