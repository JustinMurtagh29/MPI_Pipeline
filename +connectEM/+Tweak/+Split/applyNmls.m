% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
aggloFile = fullfile(rootDir, 'aggloState', 'axons_18_b.mat');
nmlDir = '/home/amotta/Desktop/test';

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

agglos = load(aggloFile);
agglos = agglos.axons(agglos.indBigAxons);

%% Find all NML files
nmlFiles = dir(fullfile(nmlDir, '*.nml'));
nmlFiles = fullfile(nmlDir, {nmlFiles(~[nmlFiles.isdir]).name});

%% Process super-agglomerates
for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    
   [~, curAggloIdName] = fileparts(curNmlFile);
    curAggloIdName = aggloIdFromFileName(curAggloIdName);
    
    curNml = slurpNml(curNmlFile);
    curAggloId = aggloIdFromNml(curNml);
    assert(curAggloId == curAggloIdName);
    
    curOut = connectEM.Tweak.applyNml(agglos(curAggloId), curNml);
    curOut = SuperAgglo.clean(curOut);
end

%% Utilities
function id = aggloIdFromFileName(fileName)
    id = regexpi(fileName, '^\d+_\w+-(?<id>\d+)', 'names', 'once');
    id = str2double(id.id);
end

function id = aggloIdFromNml(nml)
    id = regexpi( ...
        nml.parameters.experiment.description, ...
        'Agglomerate\s+(?<id>\d+)$', 'names', 'once');
    id = str2double(id.id);
end
