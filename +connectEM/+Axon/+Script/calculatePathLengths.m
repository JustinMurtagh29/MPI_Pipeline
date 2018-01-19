% Calculates the (minimum spanning tree-based) path length for each axon
% super-agglomerates.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_16_b.mat');


[outDir, outFile] = fileparts(axonFile);
outFile = sprintf('%s_meta.mat', outFile);
outFile = fullfile(outDir, outFile);
clear outDir;

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
append = {{}, {'-append'}};
append = append{1 + (exist(outFile, 'file') ~= 0)};

out = struct;
out.pathLen = pathLen;
out.pathLenInfo = info;

% write result
save(outFile, '-struct', 'out', append{:});

%% plotting
% load(outFile, 'pathLen');
pathLenUm = pathLen ./ 1E3;

figure;
histogram(pathLenUm, linspace(0, 500, 100));
title('Axonal path length distribution');
xlabel('Path length (MST) [µm]');
xlim([0, 500]);

fig = figure;
ax = axes(fig);
histogram(ax, pathLenUm, linspace(0, 500, 100));
title(ax, 'Axonal path length distribution');
xlabel(ax, 'Path length (MST) [µm]');
xlim(ax, [0, 500]);
ax.YScale = 'log';
gca.YScale = 'log';