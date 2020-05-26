% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
somaFile = fullfile(rootDir, 'aggloState', 'somata_07.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');
outFile  = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_autoSpines_v1.mat');

% Set for debugging
nmlFile = '';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

soma = load(somaFile);
dend = load(dendFile);

%% Remove somata from dendrite super-agglomerates
somata = reshape(soma.somata, [], 1);
somaAgglos = Agglo.fromSuperAgglo(somata);

out = struct;
out.dendrites = SuperAgglo.removeSegIds( ...
    param, dend.dendrites, cell2mat(somaAgglos));

%% Find whole cell for each soma
wcIds = dend.idxWholeCells(dend.indWholeCells);
wcAgglos = Agglo.fromSuperAgglo(dend.dendrites(dend.indWholeCells));

maxSegId = Seg.Global.getMaxSegId(param);
wcLUT = Agglo.buildLUT(maxSegId, wcAgglos);
somaWcIds = cellfun(@(ids) mode(nonzeros(wcLUT(ids))), somaAgglos);

% Sanity checks
% NOTE(amotta): Somata without whole cell will be NaN. Because any test for
% equality with NaN will be negative the assertion below will work even in
% this degenerate case.
assert(numel(somaWcIds) == numel(unique(somaWcIds)));
assert(all(ismember(1:numel(wcAgglos), somaWcIds)));

% Set whole cell ID to zero
somaWcIds(isnan(somaWcIds)) = 0;

%% Copying over meta data
out.parentIds = reshape(1:numel(out.dendrites), [], 1);

% Copy over indices
indFields = sort(fieldnames(dend));
indFields = indFields( ...
    startsWith(indFields, 'ind') ...
  | startsWith(indFields, 'idx'));

for curIdx = 1:numel(indFields)
    curName = indFields{curIdx};
    out.(curName) = dend.(curName);
    clear curName curVal;
end

% Remove empty super-agglomerates
mask = ~arrayfun(@(d) isempty(d.nodes), out.dendrites);
out = structfun(@(v) v(mask), out, 'UniformOutput', false);

%% Add somata
% AIS
out.indAIS = [out.indAIS; false(size(somata))];
out.idxAIS = [out.idxAIS; zeros(size(somata))];

% Whole cells
out.indWholeCells = [out.indWholeCells; false(size(somata))];
out.idxWholeCells = [out.idxWholeCells; zeros(size(somata))];

% Somata
out.indSomata = [false(size(out.dendrites)); true(size(somata))];
out.idxSomata = [zeros(size(out.dendrites)); somaWcIds];

% Big
out.indBigDends = [out.indBigDends; true(size(somata))];

% Add somata
out.dendrites = [out.dendrites; somata];
out.parentIds = [out.parentIds; zeros(size(somata))];

% Sanity checks
sizes = structfun(@size, out, 'UniformOutput', false);
sizes = cell2mat(struct2cell(sizes));
assert(size(unique(sizes, 'rows'), 1) == 1);

% TODO(amotta): Enable check when possible
out.dendrites = SuperAgglo.clean(out.dendrites, false);
out.dendAgglos = Agglo.fromSuperAgglo(out.dendrites);

%% Saving result
out.info = info;

Util.saveStruct(outFile, out);
Util.protect(outFile);

%% Debugging facilities
if ~isempty(nmlFile)
    skel = skeleton();
    skel = Skeleton.setParams4Pipeline(skel, param);
    skel = skel.setDescription(sprintf( ...
        '%s (%s)', info.filename, info.git_repos{1}.hash));

    nanSomaIds = find(isnan(somaWcIds));
    for curId = reshape(nanSomaIds, 1, [])
        skel = Skeleton.fromMST( ...
            soma.somata(curId).nodes(:, 1:3), ...
            param.raw.voxelSize, skel);
        skel.names{end} = sprintf('Soma %d', nanSomaIds);
    end

    nanWcIds = setdiff(dend.idxWholeCells, 0);
    nanWcIds(ismember(nanWcIds, somaWcIds)) = [];
    for curId = reshape(nanWcIds, 1, [])
        curDendId = find(dend.idxWholeCells == curId);

        skel = skel.addTree( ...
            sprintf('Whole cell %d', curId), ...
            dend.dendrites(curDendId).nodes(:, 1:3), ...
            dend.dendrites(curDendId).edges);
    end

    skel.write(nmlFile);
end
