% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
resultsFile = fullfile(rootDir, 'chiasmataSplitting/20171009T193744-kmb-on-axons-6c/outputs/20171030T181930_results.mat');
outputDir = '/home/amotta/Desktop/chiasma-queries';

% axons to visualize
axonIds = [321, 23565, 2206, 4566, 1388, 769, 2274];

info = Util.runInfo();

%% loading data
param = load(fullfile(rootDir, 'allParameter.mat'), 'p');
param = param.p;

data = load(resultsFile);

%%
bigAxons = data.oldAxons.indBigAxons;
bigAxons = data.oldAxons.axons(bigAxons);
bigAxons = bigAxons(axonIds);

skels = skeleton();
skels = skels.setDescription(sprintf( ...
    '%s (%s)', mfilename, info.git_repos{1}.hash));
skels = Skeleton.setParams4Pipeline(skels, param);
skels = Superagglos.buildAggloAndFlightSkels(bigAxons, skels);

for curIdx = 1:numel(skels)
    curSkelFile = sprintf(...
        '%d_axon-%d.nml', curIdx, axonIds(curIdx));
    skels(curIdx).write(fullfile(outputDir, curSkelFile));
end