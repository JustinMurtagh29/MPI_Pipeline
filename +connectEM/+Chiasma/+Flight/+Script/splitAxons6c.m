% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
axonFile = fullfile(rootDir, 'aggloState', 'axons_06_c.mat');

curDir = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171009T193744-kmb-on-axons-6c');
outputDir = fullfile(curDir, 'outputs');

% List with input files **in decreasing order of dominance**.
% Add new requery rounds to the top of this list.
dataFiles = { ...
    'requeries/20171027T102403_input-data.mat';
    'requeries/20171023T102000_input-data.mat';
    '20171018T205736_input-data.mat'};
dataFiles = fullfile(curDir, dataFiles);
clear curDir;

% run info (with above variables)
info = Util.runInfo();

%% loading parameter file
load(fullfile(rootDir, 'allParameter.mat'), 'p');

%% splitting chiasmata

%%
outFile = sprintf('%s_results.mat', datestr(now, 30));
outFile = fullfile(outputDir, outFile);
Util.saveStruct(outFile, out);