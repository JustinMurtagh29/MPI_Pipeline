% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outFile = '/tmpscratch/amotta/l4/2018-05-08-ending-query-evaluation/axon-state_60.mat';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Find nml Files
skelDirs = connectEM.setQueryState('6.0');
skelDirs = fullfile('/gaba', reshape(skelDirs, [], 1));

nmlFiles = cellfun( ...
    @(skelDir) dir(fullfile(skelDir, '*.nml')), ...
    skelDirs, 'UniformOutput', false);
nmlFiles = arrayfun( ...
    @(f) fullfile(f.folder, f.name), ...
    cat(1, nmlFiles{:}), 'UniformOutput', false);

%% Processing
nmlT = table;
nmlT.filePath = nmlFiles;
nmlT.pathLenNm = nan(height(nmlT), 1);
nmlT.timeMs = nan(height(nmlT), 1);

tic;
for curIdx = 1:height(nmlT)
    curFilePath = nmlT.filePath{curIdx};
    
    try
       [nmlT.pathLenNm(curIdx), nmlT.timeMs(curIdx)] = ...
            forNmlFile(curFilePath, param.raw.voxelSize);
    catch
        warning('Error in %s', curFilePath);
    end
    
    Util.progressBar(curIdx, height(nmlT));
end

%% Save result
out = struct;
out.nmlT = nmlT;
out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);

%% Utility
function [lenNm, timeMs] = forNmlFile(path, voxelSize)
    nml = slurpNml(path);
    nodes = NML.buildNodeTable(nml);
    nodes.coord = nodes.coord .* voxelSize;
    
    timeMs = sort(nodes.time);
    timeMs = timeMs(end) - timeMs(2);
    
    lenNm = cell2mat(cellfun( ...
        @(e) horzcat(e.source, e.target), ...
        nml.things.edges, 'UniformOutput', false));
   [~, lenNm] = ismember(lenNm, nodes.id);
   
    lenNm = ...
        nodes.coord(lenNm(:, 1)) ...
      - nodes.coord(lenNm(:, 2));
    lenNm = sum(lenNm .* lenNm, 2);
    lenNm = sum(sqrt(lenNm));
end
