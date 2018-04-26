% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-classified_spine-syn-clust.mat');

% This corresponds to the axon shown in panel 3c. This was axon 630 in
% commectome `connectome_axons_18_a_ax_spine_syn_clust`.
axonId = 628;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, ~, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[classConnectome, targetClasses] = ...
    connectEM.Connectome.buildClassConnectome(conn);
