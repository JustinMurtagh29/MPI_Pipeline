% Calculates the (minimum spanning tree-based) path length for each axon
% super-agglomerates.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');

info = Util.runInfo();

%% loading data
fprintf('Loading data... ');
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

agglos = load(axonFile, 'axons', 'indBigAxons');
agglos = agglos.axons(agglos.indBigAxons);
fprintf('done!\n');

%% calculate path length
fprintf('Calculating path lengths... ');
pathLen = Superagglos.mstLength(agglos, param.raw.voxelSize);
fprintf('done!\n');

%% write output
[outDir, outFile] = fileparts(axonFile);
outFile = sprintf('%s_meta.mat', outFile);
outFile = fullfile(outDir, outFile);

append = {{'-append'}, {}};
append = append{1 + logical(exist(outFile, 'file'))};

out = struct;
out.pathLen = pathLen;
out.pathLenInfo = info;

% write result
save(outFile, '-struct', 'out', append{:});