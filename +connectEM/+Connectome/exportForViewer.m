% This code is based on
%   +connectEM/+Figure/connectome.m
%   a8e72777a9b77d36addbe7f3f701ac014192585f
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
outDir = '/tmpscratch/amotta/l4/2019-08-23-export-to-connectome-viewer';

axonClasses = { ...
    'Corticocortical', ...
    'Thalamocortical', ...
    'Other', ...
    'Inhibitory'};
dendClasses = { ...
    'Soma', ...
    'ProximalDendrite', ...
    'SmoothDendrite', ...
    'ApicalDendrite', ...
    'AxonInitialSegment', ...
    'OtherDendrite'};

minSynCount = 10;
runId = datestr(now, 30)

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[rawConn, rawSyn] = connectEM.Connectome.load(param, connFile);

rawAxons = load(rawConn.info.param.axonFile);
rawAxons = rawAxons.axons;

rawDends = load(rawConn.info.param.dendriteFile);
rawDends = rawDends.dendrites;

%% Filter neurites
conn = rawConn;

% Only plot neurites with at least `minSynCount` synapses. Notable
% exception are AIS, which we plot irrespectively of their synapse number.
axonMeta = conn.axonMeta;
axonMeta(axonMeta.synCount < minSynCount, :) = [];

dendMeta = conn.denMeta;
dendMeta( ...
    dendMeta.targetClass ~= 'AxonInitialSegment' ...
  & dendMeta.synCount < minSynCount, :) = [];

% Rename "Somata" to "Soma"
somaMask = dendMeta.targetClass == 'Somata';
dendMeta.targetClass(somaMask) = 'Soma';

% Rename "WholeCell" to "ProximalDendrite"
proxDendMask = dendMeta.targetClass == 'WholeCell';
dendMeta.targetClass(proxDendMask) = 'ProximalDendrite';

conn = conn.connectome;
conn(~ismember(conn.edges(:, 1), axonMeta.id), :) = [];
conn(~ismember(conn.edges(:, 2), dendMeta.id), :) = [];

%% Group by classes
rng(0);
curIds = randperm(numel(axonMeta.id));
axonMeta = sortrows(axonMeta, 'id');
axonMeta = axonMeta(curIds, :);

[~, axonMeta.classId] = ismember( ...
    axonMeta.axonClass, axonClasses);

assert(all(axonMeta.classId));
axonMeta = sortrows(axonMeta, 'classId');

rng(0);
curIds = randperm(numel(dendMeta.id));
dendMeta = sortrows(dendMeta, 'id');
dendMeta = dendMeta(curIds, :);

[~, dendMeta.classId] = ismember( ...
    dendMeta.targetClass, dendClasses);

assert(all(dendMeta.classId));
dendMeta = sortrows(dendMeta, 'classId');

%% Restrict super-agglomerates
clear cur*;

curCleanSuperAgglos = @(s) rmfield( ...
    s, setdiff(fieldnames(s), {'edges', 'nodes'}));
axons = curCleanSuperAgglos(rawAxons);
dends = curCleanSuperAgglos(rawDends);

[curAxonIds, ~, axonMeta.skelId] = unique(axonMeta.parentId);
axons = axons(curAxonIds(:));

[curDendIds, ~, dendMeta.skelId] = unique(dendMeta.parentId);
dends = dends(curDendIds(:));

%% Build output
clear cur*;

out = struct;
out.axonClasses = axonClasses(:);
curAxonVars = {'id', 'skelId', 'classId'};
out.axonMeta = axonMeta(:, curAxonVars);
out.axonSkels = axons(:);

out.dendClasses = dendClasses(:);
curDendVars = {'id', 'skelId', 'classId', 'cellId', 'isInterneuron'};
out.dendMeta = dendMeta(:, curDendVars);
out.dendSkels = dends(:);

out.conn = conn;
out.conn.Properties.VariableNames = {'edge', 'synIds'};

out.info = info;

% Simplify format
curFields = fieldnames(out);
for curIdx = 1:numel(curFields)
    curData = out.(curFields{curIdx});
    
    if istable(curData)
        % Convert table to structure of arrays
        curData = table2struct(curData, 'ToScalar', true);
    elseif isstruct(curData) && numel(curData) > 1
        % Convert array of structures to structure of array
        curSubFields = fieldnames(curData);
        
        curData = cellfun( ...
            @(n) reshape({curData.(n)}, [], 1), ...
            curSubFields, 'UniformOutput', false);
        curData = cell2struct(curData, curSubFields);
    end
    
    out.(curFields{curIdx}) = curData;
end

%% Save result
clear cur*;

curOutName = sprintf('%s_bundle.mat', runId);
curOutPath = fullfile(outDir, curOutName);

save(curOutPath, '-struct', 'out', '-v6');
Util.protect(curOutPath);
