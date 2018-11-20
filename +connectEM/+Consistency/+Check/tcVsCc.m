% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
connLinFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn] = connectEM.Connectome.load(param, connFile);

% HACK(amotta): The connectome used `axons_19_a_partiallySplit_v2.mat`, but
% the parent IDs were only added later in version 4...
axonMeta = load(strrep(conn.info.param.axonFile, '_v2.mat', '_v4.mat'));
conn.axonMeta.unsplitParentId = axonMeta.parentIds(conn.axonMeta.parentId);
clear axonMeta;

[connLin, synLin] = connectEM.Connectome.load(param, connLinFile);

axonMeta = load(connLin.info.param.axonFile);
connLin.axonMeta.unsplitParentId = ...
    axonMeta.parentIds(connLin.axonMeta.parentId);
clear axonMeta;

%% Translate axon classes from partially split to fully linearized version
% NOTE(amotta): The partially split and the linearized version of the axons
% were both derived from axons 19a. To translate the axon classes from the
% partially split to the fully linearized axons, we thus have to do a
% detour via the unsplit version.
%   Because there's no bijective mapping between partially split and fully
% linearized axons, we only look at the axons 19a which were uniquely
% classified in the partially split version.
clear cur*;

% Find parents of axons classified as TC / CC
curTcParentIds = conn.axonMeta.unsplitParentId( ...
    conn.axonMeta.axonClass == 'Thalamocortical');
curCcParentIds = conn.axonMeta.unsplitParentId( ...
    conn.axonMeta.axonClass == 'Corticocortical');

% Check if there exist siblings of TC axons that were classified
% differently (i.e., as corticocortical or inhibitory; other is okay). If
% yes, we do not trust or use these TC axons.
curTcSiblingsOkay = accumarray( ...
    conn.axonMeta.unsplitParentId, ...
    double(conn.axonMeta.axonClass), [], ...
    @(types) isequal(setdiff(types, 4), 2)); % 2 = TC, 4 = other
curTcSiblingsOkay = curTcSiblingsOkay(curTcParentIds);
curTcParentIds = curTcParentIds(curTcSiblingsOkay);

curCcSiblingsOkay = accumarray( ...
    conn.axonMeta.unsplitParentId, ...
    double(conn.axonMeta.axonClass), [], ...
    @(types) isequal(setdiff(types, 4), 1)); % 1 = CC, 4 = other
curCcSiblingsOkay = curCcSiblingsOkay(curCcParentIds);
curCcParentIds = curCcParentIds(curCcSiblingsOkay);

% Sanity check
assert(isempty(intersect( ...
    curTcParentIds, curCcParentIds)));

% Assign axon classes to linearized axons
connLin.axonMeta.axonClass(:) = 'Other';
connLin.axonMeta.axonClass(ismember( ...
    connLin.axonMeta.unsplitParentId, ...
    curTcParentIds)) = 'Thalamocortical';
connLin.axonMeta.axonClass(ismember( ...
    connLin.axonMeta.unsplitParentId, ...
    curCcParentIds)) = 'Corticocortical';
