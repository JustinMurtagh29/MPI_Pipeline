% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuratiom
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');

pairNames = { ...
    % farther than 20 µm
    '02_axon-16935_dendrite_11105';
    '03_axon-28035_dendrite_11281';
    '04_axon-52046_dendrite_450';
    '05_axon-63838_dendrite_11252';
    '07_axon-38403_dendrite_1572';
    '08_axon-43563_dendrite_11226';
    '09_axon-4964_dendrite_11114';
    '13_axon-4270_dendrite_5796';
    '14_axon-40081_dendrite_11257';
    '15_axon-48021_dendrite_4468';
    '16_axon-28048_dendrite_5315';
    '17_axon-5057_dendrite_2836';
    '18_axon-25201_dendrite_1093';
    '20_axon-16318_dendrite_2417';
    '23_axon-50921_dendrite_5179';
    '24_axon-67244_dendrite_11270';
    % closer than 5 µm
    '01_axon-55107_dendrite_272';
    '02_axon-65041_dendrite_962';
    '04_axon-8046_dendrite_11245';
    '07_axon-64124_dendrite_2589';
    '09_axon-20537_dendrite_5060';
    '11_axon-67527_dendrite_11260';
    '12_axon-32601_dendrite_3841';
    '13_axon-23876_dendrite_2644';
    '14_axon-45881_dendrite_8835';
    '16_axon-50397_dendrite_1939'};

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

[~, synToSynFile] = fileparts(connFile);
synToSynFile = sprintf('%s_synToSynDists.mat', synToSynFile);
synToSynFile = fullfile(fileparts(connFile), synToSynFile);
synToSyn = load(synToSynFile);

%% Prepare neurite pairs
pairPattern = '\d+_axon-(?<preAggloId>\d+)_dendrite_(?<postAggloId>\d+)';

pairConfig = regexp(pairNames, pairPattern, 'names');
pairConfig = cellfun(@(pair) structfun( ...
    @str2num, pair, 'UniformOutput', false), pairConfig);
pairConfig = struct2table(pairConfig);

% Find corresponding synapse pairs
pairSynT = synT;
pairSynT.id = reshape(1:height(pairSynT), [], 1);

pairSynT(~pairSynT.isSpine, :) = [];

[~, pairSynT.pairId] = ismember(pairSynT(:, { ...
    'preAggloId', 'postAggloId'}), pairConfig, 'rows');
pairSynT(~pairSynT.pairId, :) = [];

pairSynT = sortrows(pairSynT, 'pairId');
pairConfig.synIdPairs = transpose(reshape(pairSynT.id, 2, []));

pairConfig = table2struct(pairConfig, 'ToScalar', true);
pairConfig.title = sprintf( ...
    'Proofread synapse pairs (n = %d)', ...
    size(pairConfig.synIdPairs, 1));

%% Plot
connectEM.Consistency.plotVariabilityVsDistance( ...
    synT, synToSyn, pairConfig, 'info', info);
