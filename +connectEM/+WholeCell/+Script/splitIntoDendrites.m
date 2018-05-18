% Because of a damn webKNOSSOS bug we now have to manually reconstruct
% which NML files corresponds to which whole cell...
% 
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto-and-manual.mat');
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

[~, wcSplitFile] = pdist2( ...
    param.raw.voxelSize .* nmlFirstNodes, ...
    param.raw.voxelSize .* aggloFirstNodes, ...
    'squaredeuclidean', 'Smallest', 1);

% Make sure this crap is bijective
assert(numel(wcSplitFile) == numel(unique(wcSplitFile)));
wcSplitFile = nmlFiles(wcSplitFile);
