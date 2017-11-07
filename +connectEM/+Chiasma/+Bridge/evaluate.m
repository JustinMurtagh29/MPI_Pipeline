% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
chiasmataFile = fullfile( ...
    rootDir, 'chiasmataSplitting', ...
    '20171104T181213-on-axons-10a-plus-10kE3a', ...
    '20171104T184018_chiasmata.mat');

outputDir = '/home/amotta/Desktop/random-bridges';
debug = true;

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

% shuffle order
chiasmaCount = size(chiasma, 1);

%%
paramForChiasma = transpose(horzcat( ...
    fieldnames(chiasmaParam), struct2cell(chiasmaParam)));
paramForChiasma = Util.modifyStruct(param, paramForChiasma{:});

results = table;
results.nrPartitions = nan(chiasmaCount, 1);
results.nrBridges = nan(chiasmaCount, 1);

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

tic
for curIdx = 1:chiasmaCount
    curAxonId = chiasma.axonId(curIdx);
    curAxon = axons(curAxonId);
    
    curChiasmata = chiasmata{curAxonId};
    curChiasmaId = chiasma.chiasmaId(curIdx);
    
    try
       [isBridge, edges] = connectEM.Chiasma.Bridge.findBridges( ...
            param, chiasmaParam, curAxon, curChiasmata, curChiasmaId);

        results.nrPartitions(curIdx) = sum(any(isBridge, 1));
        results.nrBridges(curIdx) = sum(isBridge(:));
    catch
        warning('FIX ME!');
        continue;
    end
    
    % run the rest only in debug mode
    if ~exist('debug', 'var') || ~debug; continue; end
    
    % restrict debugging to most interesting cases
    if results.nrPartitions(curIdx) ~= 1; continue; end
    
    curComments = repmat({''}, size(curAxon.nodes, 1), 1);
    curComments{curChiasmata.ccCenterIdx(curChiasmaId)} = 'Chiasma';
    curComments(curChiasmata.queryIdx{curChiasmaId}) = {'Chiasma exit'};
    
    curSkel = skeleton();
    curSkel = curSkel.addTree( ...
        sprintf('Axon %d', curAxonId), curAxon.nodes(:, 1:3), ...
        curAxon.edges, [], [], curComments);
    
    curBridgeColor = [0, 1, 0, 1];
    curBridges = edges(any(isBridge, 2), :);
    
    for curBridgeIdx = 1:size(curBridges, 1)
        curNodes = curBridges(curBridgeIdx, :);
        curNodes = curAxon.nodes(curNodes, 1:3);
        
        curSkel = curSkel.addTree( ...
            sprintf('Bridge %d', curBridgeIdx), ...
            curNodes, [1, 2], curBridgeColor);
    end
    
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    curSkel.write(fullfile(outputDir, sprintf( ...
        'axon-%d-chiasma-%d.nml', curAxonId, curChiasmaId)));
end
toc

%% show results
fprintf('\n');
fprintf('# chiasmata: %d\n', chiasmaCount);
fprintf('# splittable chiasmata: %d\n', sum(results.nrPartitions > 0));
fprintf('# chiasmata with unique split: %d\n', sum(results.nrPartitions == 1));