% TODO(amotta):
% [ ] Should we only export axons with ≥ 10 synapses? What about dendrites?
% [ ] All synapses or only the ones collected by agglomerates?
% [ ] Synapse positions
% [ ] Axon spine interfaces
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified__20190308T100942_asiT.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
outDir = '/tmpscratch/amotta/l4/2019-10-09-axons-dendrites-connectome-in-hdf5';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);

axons = load(conn.info.param.axonFile);
dendrites = load(conn.info.param.dendriteFile);

shAgglos = load(shFile);
shAgglos = shAgglos.shAgglos;

asiT = load(asiFile);
asiT = asiT.asiT;

%% Preprocess data
clear cur*;

% Rename "Somata" to "Soma"
curMask = conn.denMeta.targetClass == 'Somata';
conn.denMeta.targetClass(curMask) = 'Soma';

% Rename "WholeCell" to "ProximalDendrite" or "SmoothDendrite"
curInMask = conn.denMeta.isInterneuron;
curWcMask = conn.denMeta.targetClass == 'WholeCell';
conn.denMeta.targetClass(curWcMask &  curInMask) = 'SmoothDendrite';
conn.denMeta.targetClass(curWcMask & ~curInMask) = 'ProximalDendrite';

axons = axons.axons(conn.axonMeta.parentId);
dendrites = dendrites.dendrites(conn.denMeta.parentId);

%% Build info string
clear cur*;
curIsDirtyStrs = {'', ' (dirty)'};

infoStr = cell(0, 1);
infoStr{end + 1} = sprintf('%s', info.filename);

for curRepoId = 1:numel(info.git_repos)
    curRepo = info.git_repos{curRepoId};
    infoStr{end + 1} = sprintf( ...
        '%s %s%s', curRepo.remote, curRepo.hash, ...
        curIsDirtyStrs{2 - isempty(curRepo.diff)}); %#ok
end

infoStr{end + 1} = sprintf( ...
    '%s@%s. MATLAB %s. %s', ...
    info.user, info.hostname, ...
    info.matlab_version, info.time);
infoStr = strjoin(infoStr, newline);

%% Export axon and dendrite reconstructions
clear cur*;

% Axons
curOutFile = fullfile(outDir, 'axons.hdf5');
categoricalToHdf5(curOutFile, '/axonClasses', conn.axonMeta.axonClass);
agglosToHdf5(curOutFile, '/axonAgglomerates', conn.axons);
superAgglosToHdf5(curOutFile, '/axonSkeletons', axons);
h5writeatt(curOutFile, '/', 'info', infoStr);
Util.protect(curOutFile);

% Dendrites
curOutFile = fullfile(outDir, 'dendrites.hdf5');
categoricalToHdf5(curOutFile, '/targetClasses', conn.denMeta.targetClass);
agglosToHdf5(curOutFile, '/dendriteAgglomerates', conn.dendrites);
superAgglosToHdf5(curOutFile, '/dendriteSkeletons', dendrites);
agglosToHdf5(curOutFile, '/spineHeadAgglomerates', shAgglos);
h5writeatt(curOutFile, '/', 'info', infoStr);
Util.protect(curOutFile);

%% Export connectome
clear cur*;
%{
curSyn = syn.synapses(:, {'presynId', 'postsynId'});
curSyn.Properties.VariableNames = {'preSegIds', 'postSegIds'};
curSyn = table2struct(curSyn, 'ToScalar', false);

curRepIds = cellfun(@numel, conn.connectome.synIdx);
curRepIds = repelem(reshape(1:numel(curRepIds), [], 1), curRepIds(:));

curConn = table;
curConn.synId = cell2mat(conn.connectome.synIdx);
curConn.preAxonId = conn.connectome.edges(curRepIds, 1);
curConn.postDendId = conn.connectome.edges(curRepIds, 2);
curConn = table2struct(curConn, 'ToScalar', true);

curOutFile = fullfile(outDir, 'connectome.hdf5');

categoricalToHdf5(curOutFile, '/synapseTypes', syn.synapses.type);
structToHdf5(curOutFile, '/synapses', curSyn);
%}

%% Utilities
function superAgglosToHdf5(outFile, group, agglos)
    agglos = rmfield(agglos, setdiff( ...
        fieldnames(agglos), {'nodes', 'edges'}));
    
    for curIdx = 1:numel(agglos)
        agglos(curIdx).nodes = reshape( ...
            agglos(curIdx).nodes(:, 1:3), [], 3);
        agglos(curIdx).edges = reshape(unique( ...
            sort(agglos(curIdx).edges, 2), 'rows'), [], 2);
    end
    
    structToHdf5(outFile, group, agglos, true);
end

function agglosToHdf5(outFile, group, agglos)
    assert(all(cellfun(@(ids) all(isfinite(ids)), agglos)));
    assert(all(cellfun(@(ids) all(ids <= intmax('uint32')), agglos)));
    
    agglos = cellfun( ...
        @(ids) uint32(ids(:)), ...
        agglos(:), 'UniformOutput', false);
    cellToHdf5(outFile, group, agglos);
end

function cellToHdf5(file, group, cellArray)
    assert(isvector(cellArray));
    for curIdx = 1:numel(cellArray)
        curData = cellArray{curIdx};
        curDset = fullfile(group, num2str(curIdx));
        numericToHdf5(file, curDset, curData);
    end
end

function structToHdf5(file, group, struct, forceArray)
    if ~exist('forceArray', 'var') || isempty(forceArray)
        forceArray = false;
    end
    
    if ~isscalar(struct) || forceArray
        for curIdx = 1:numel(struct)
            curGroup = fullfile(group, num2str(curIdx));
            structToHdf5(file, curGroup, struct(curIdx));
        end
    else
        names = fieldnames(struct);
        values = struct2cell(struct);
        assert(isequal(numel(names), numel(values)));

        for curIdx = 1:numel(names)
            curDset = fullfile(group, names{curIdx});
            numericToHdf5(file, curDset, values{curIdx});
        end
    end
end

function categoricalToHdf5(outFile, dset, cats)
    assert(iscategorical(cats));
    catNames = categories(cats);
    
    assert(numel(catNames) <= intmax('uint8'));
    numericToHdf5(outFile, dset, uint8(cats));
    
    for curIdx = 1:numel(catNames)
        curName = lower(catNames{curIdx});
        h5writeatt(outFile, dset, curName, uint8(curIdx));
    end
end

function numericToHdf5(file, dset, data)
    assert(isnumeric(data));
    sz = size(data);
    
    % NOTE(amotta): Remove trailing singleton dimension
    if numel(sz) == 2 && sz(2) == 1; sz = sz(1); end
    
    h5create( ...
        file, dset, sz, ...
        'Datatype', class(data), ...
        'Chunksize', sz, ...
        'Deflate', 9);
    h5write(file, dset, data);
end
