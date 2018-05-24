% Because of a damn webKNOSSOS bug we now have to manually reconstruct
% which NML files corresponds to which whole cell...
% 
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendFile = fullfile(rootDir, 'aggloState', ...
    'dendrites_wholeCells_02_v3_auto-and-manual.mat');
nmlDir = fullfile(fileparts(mfilename('fullpath')), 'annotations');

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

dend = load(dendFile);

%% Loading NML files
nmlFiles = dir(fullfile(nmlDir, '*.nml'));
nmlFiles = fullfile(nmlDir, {nmlFiles.name});

[wcIds, aggloIds] = setdiff(dend.idxWholeCells, 0);
agglos = dend.dendrites(aggloIds);

aggloFirstNodes = cell2mat(arrayfun( ...
    @(agglo) agglo.nodes(1, 1:3), ...
    agglos, 'UniformOutput', false));

nmlFirstNodes = nan(numel(nmlFiles), 3);
for curFileIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curFileIdx};
    curNml = slurpNml(curNmlFile);
    
    curNodes = NML.buildNodeTable(curNml);
    curNodes.coord = curNodes.coord + 1;
    
    nmlFirstNodes(curFileIdx, :) = ...
        curNodes.coord(curNodes.id == 1, :);
end

[~, nmlAggloIds] = pdist2( ...
    param.raw.voxelSize .* aggloFirstNodes, ...
    param.raw.voxelSize .* nmlFirstNodes, ...
    'squaredeuclidean', 'Smallest', 1);
nmlAggloIds = aggloIds(nmlAggloIds);

% Make sure this crap is bijective
assert(numel(nmlAggloIds) ...
    == numel(unique(nmlAggloIds)));

%% Patch NML files
for curIdx = 1:numel(nmlFiles)
    curNmlFile = nmlFiles{curIdx};
    curAggloId = nmlAggloIds(curIdx);
    
    curDesc = sprintf('Agglomerate %d', curAggloId);
    curRegexp = sprintf('s/description=".*"/description="%s"/g', curDesc);
    curCommand = sprintf('sed -i ''%s'' "%s"', curRegexp, curNmlFile);
    curError = system(curCommand);
    assert(~curError);
end
