% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-04_dendrites-wholeCells-autoSpines-v1-classified-v2_SynapseAgglos-autoPreRobo-v1-classified.mat');
splitAxonFile = fullfile(rootDir, 'aggloState', 'axons_04_linearized.mat');
outFile = fullfile(rootDir, 'aggloState', 'axons_04_partiallySplit_v1.mat');

minSynPre = 10;

% Set file path to generate debug NML file
debugNmlFile = '/gaba/u/amotta/problematic-axons-04_v1.nml';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

% Loading axons
axons = load(conn.info.param.axonFile);

% NOTE(amotta): This early in the axon reconstruction, we haven't set the
% `endings` and `solvedChiasma` fields yet. So, let's just set them to
% correctly shaped null values.
for curIdx = 1:numel(axons.axons)
    axons.axons(curIdx).endings = zeros(0, 1);
    axons.axons(curIdx).solvedChiasma = ...
        false(size(axons.axons(curIdx).nodes, 1), 1);
end

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
    skel = Skeleton.setDescriptionFromRunInfo(skel, info);
    
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
