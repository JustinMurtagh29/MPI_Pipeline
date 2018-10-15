% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/amotta/l23/2018-10-09-mrnet-pipeline-run';
outDir = '/tmpscratch/amotta/l23/2018-10-13-vessel-and-nuclei-mask';

mag = [8, 8, 4];

nmlFile = fileparts(mfilename('fullpath'));
nmlFile = fullfile(nmlFile, '+Data', 'nuclei_merger-mode.nml');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading parameters
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

%% Loading segmentation
segParam = param.seg;
segParam.root = fullfile( ...
    fileparts(segParam.root(1:(end - 1))), ...
    sprintf('%d-%d-%d', mag));

magBox = ceil(param.bbox ./ mag(:));
seg = loadSegDataGlobal(segParam, magBox);

%% Loading merger mode tracings
nml = slurpNml(nmlFile);
nodes = NML.buildNodeTable(nml);

nodes.coord = nodes.coord + 1;
nodes.coord = nodes.coord - param.bbox(:, 1)' + 1;
nodes.coord = ceil(nodes.coord ./ mag);

assert(all(all(nodes.coord >= 1)));
assert(all(all(nodes.coord <= size(seg))));

nodes.segId = Util.sub2ind(size(seg), nodes.coord);
nodes.segId = seg(nodes.segId);

%% Build mask from segments
maxSegId = Seg.Global.getMaxSegId(param);

lut = false(1 + maxSegId, 1);
lut(1 + nodes.segId) = true;

seg = seg + 1;
nucMask = lut(seg);
clear seg;

%% Fill holes
nucMask = padarray(nucMask, [1, 1, 1], true);

nucMask = imclose(nucMask, strel('sphere', 1));
nucMask(1:2, 1:2, round(size(nucMask, 3) / 2)) = false;
nucMask = imfill(nucMask, 'holes');

nucMask([1, end], :, :) = [];
nucMask(:, [1, end], :) = [];
nucMask(:, :, [1, end]) = [];

%% Save result
outFile = fullfile(outDir, 'nuclei_v1.mat');
Util.save(outFile, nucMask, info);
