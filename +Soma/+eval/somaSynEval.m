% evaluation of soma synapses
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo();
outFolder = ['/gaba/u/mberning/results/pipeline/20170217_ROI/' ...
    'soma_BS/SomaSynEval/'];
if ~exist(outFolder, 'dir')
    mkdir(outFolder)
end


%% select three random center somata for synapse annotation by hiwis

p = Gaba.getSegParameters('ex145_ROI2017');
m = load(p.agglo.somaFile);
somaAgglos = m.somaAgglos(m.centerSomaIdx, 1);

% get 3 random soma
n = 3;
rng('shuffle')
idx = randperm(length(somaAgglos));
idx = idx(1:3)';

% get center of mass
m = load(p.svg.segmentMetaFile, 'point');
point = m.point';
somaSeeds = cell2mat(cellfun(@(x)round(mean(point(x, :), 1)), ...
    somaAgglos(idx), 'uni', 0));

% results
seeds = table(idx, somaSeeds, 'VariableNames', {'SomaId', 'Seed'});

outFile = fullfile(outFolder, 'SomaSynEval_seeds.mat');
if ~exist(outFile, 'file')
    Util.log('Saving seeds to %s.', outFile);
    save(outFile, 'seeds', 'info');
else
    Util.log('Output file %s already exists.');
end

