% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
chiasmataFile = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a', ...
    '20171104T184018_chiasmata.mat');

outputDir = '/home/amotta/Desktop/random-chiasmata';

%% load data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

data = load(chiasmataFile);
axonFile = data.info.param.axonFile;
chiasmaParam = data.info.param.chiasmaParam;
chiasmata = data.chiasmata;
clear data;

axons = load(axonFile);
axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% select examples
chiasma = table;
chiasma.axonId = repelem((1:numel(chiasmata))', ...
    cellfun(@(c) numel(c.ccCenterIdx), chiasmata));
chiasma.chiasmaId = cell2mat(cellfun(@(c) ...
    (1:numel(c.ccCenterIdx))', chiasmata, 'UniformOutput', false));

% restrict to 4-fold chiasmata
chiasma.nrExits = arrayfun(@(a, c) ...
    chiasmata{a}.nrExits(chiasmata{a}.ccCenterIdx(c)), ...
    chiasma.axonId, chiasma.chiasmaId);
chiasma(chiasma.nrExits ~= 4, :) = [];
chiasma.nrExits = [];

chiasmaCount = size(chiasma, 1);

%%
paramForChiasma = transpose(horzcat( ...
    fieldnames(chiasmaParam), struct2cell(chiasmaParam)));
paramForChiasma = Util.modifyStruct(param, paramForChiasma{:});

results = table;
results.nrPartitions = nan(chiasmaCount, 1);
results.nrBridges = nan(chiasmaCount, 1);

tic
for curIdx = 1:size(chiasma, 1)
    curAxon = axons(chiasma.axonId(curIdx));
    curChiasmata = chiasmata{chiasma.axonId(curIdx)};
    curChiasmaId = chiasma.chiasmaId(curIdx);
    
   [isBridge, edges] = connectEM.Chiasma.Bridge.findBridges( ...
        param, chiasmaParam, curAxon, curChiasmata, curChiasmaId);
    
    results.nrPartitions(curIdx) = sum(any(isBridge, 1));
    results.nrBridges(curIdx) = sum(isBridge(:));
end
toc

%% show results
fprintf('\n');
fprintf('# chiasmata: %d\n', chiasmaCount);
fprintf('# splittable chiasmata: %d\n', sum(results.nrPartitions > 0));
fprintf('# chiasmata with unique split: %d\n', sum(results.nrPartitions == 1));