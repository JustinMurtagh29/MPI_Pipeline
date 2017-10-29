% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_08_a.mat');

% Some parameter for algorithm
%{
chiParam = struct;
chiParam.sphereRadiusOuter = 10000; % in nm
chiParam.sphereRadiusInner = 2500; % in nm
chiParam.minNodeDist = 3000; % in nm
chiParam.clusterSize = 3000; % in nm
%}

% old config
chiParam = struct;
chiParam.sphereRadiusOuter = 10000; % in nm
chiParam.sphereRadiusInner = 1000; % in nm
chiParam.minNodeDist = 2000; % in nm
chiParam.clusterSize = 2000; % in nm

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);

axons = load(axonFile, 'axons', 'indBigAxons');
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% add chiasma parameters
chiParam = cat(2, fieldnames(chiParam), struct2cell(chiParam));
chiParam = reshape(transpose(chiParam), 1, []);

param = Util.modifyStruct(param, chiParam{:});

%% find chiasmata with mergers
%{
% find christians "huge merger" axons
axonLUT = Agglo.buildLUT(maxSegId, Superagglos.getSegIds(axons));
theseAxonIds = axonLUT(segIds);
%}

theseAxonIds = [1872; 3352];
theseAxons = axons(theseAxonIds);
