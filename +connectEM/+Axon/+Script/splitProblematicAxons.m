% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
splitAxonFile = fullfile(rootDir, 'aggloState', 'axons_19_a_linearized.mat');
outFile = fullfile(rootDir, 'aggloState', 'axons_19_a_partiallySplit_v3.mat');

% Set file path to generate debug NML file
debugNmlFile = '';

minSynPre = 10;

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

% Loading axons
axons = load(conn.info.param.axonFile);
splitAxons = load(splitAxonFile);

% Sanity check
chiasmata = load(splitAxons.info.param.chiasmaFile, 'info');
assert(isequal(conn.info.param.axonFile, chiasmata.info.param.axonFile));

%% Find axons to split
toSplitIds = conn.axonMeta.parentId( ...
    conn.axonMeta.synCount >= minSynPre ...
  & conn.axonMeta.fullPriSpineSynFrac > 0.2 ...
  & conn.axonMeta.fullPriSpineSynFrac < 0.5);

[~, replaceAxonByIds] = ismember( ...
    splitAxons.parentIds, toSplitIds);
replaceAxonByIds = accumarray( ...
    nonzeros(replaceAxonByIds), ...
    find(replaceAxonByIds), ...
   [numel(toSplitIds), 1], ...
    @(ids) {ids});

replacedAxonIds = repelem( ...
    toSplitIds, cellfun(@numel, replaceAxonByIds));

assert(~any(cellfun(@isempty, replaceAxonByIds)));
replaceByIds = cell2mat(replaceAxonByIds);

%% Debug
if ~isempty(debugNmlFile)
    exportIds = toSplitIds(1:min(10, numel(toSplitIds)));
    numDigits = ceil(log10(1 + numel(exportIds)));
    
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));
    
    for curIdx = 1:numel(exportIds)
        curId = exportIds(curIdx);
        
        curAgglo = axons.axons(curId);
        curSplitAgglos = splitAxons.axons(replaceAxonByIds{curIdx});
        curTreeIds = skel.numTrees() + (1:(numel(curSplitAgglos) + 1));
       
        skel = Superagglos.toSkel(curAgglo, skel);
        skel.names{end} = 'Original';
        
        skel = Superagglos.toSkel(curSplitAgglos, skel);
        skel.names(curTreeIds(2:end)) = arrayfun( ...
            @(id) sprintf('Component %d', id), ...
            1:(numel(curTreeIds(2:end))), ...
            'UniformOutput', false);
        
        curGroupName = sprintf('%0*d. Axon %d', numDigits, curIdx, curId);
       [skel, curGroupId] = skel.addGroup(curGroupName);
        skel = skel.addTreesToGroup(curTreeIds, curGroupId);
    end
    
    skel.write(debugNmlFile);
end

%% Build output
out = struct;
out.axons = axons.axons;
out.indBigAxons = axons.indBigAxons;
out.parentIdsSplit = zeros(size(out.axons));
out.parentIds = reshape(1:numel(out.axons), [], 1);

% Remove axons which we want to split
out.axons(toSplitIds) = [];
out.indBigAxons(toSplitIds) = [];
out.parentIdsSplit(toSplitIds) = [];
out.parentIds(toSplitIds) = [];

% Addend split versions
out.axons = [out.axons; splitAxons.axons(replaceByIds)];
out.indBigAxons = [out.indBigAxons; splitAxons.indBigAxons(replaceByIds)];
out.parentIdsSplit = [out.parentIdsSplit; replaceByIds];
out.parentIds = [out.parentIds; replacedAxonIds];

out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);
