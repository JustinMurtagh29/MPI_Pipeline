% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
somaFile = fullfile(rootDir, 'aggloState', 'somata_06.mat');
dendFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v2_auto-and-manual.mat');

% Set for debugging
nmlFile = '/home/amotta/Desktop/somata-debug.nml';

info = Util.runInfo();

%% Loading data
Util.log('Loading data');

param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

soma = load(somaFile);
dend = load(dendFile);

%% Remove somata from dendrite super-agglomerates
Util.log('Removing somata from dendrites');
somaAgglos = Agglo.fromSuperAgglo(soma.somata);

out = struct;
out.dendrites = SuperAgglo.removeSegIds( ...
    param, dend.dendrites, cell2mat(somaAgglos));
out.dendrites = SuperAgglo.clean(out.dendrites, false);

%% Add somata to list of postsynaptic processes
wcIds = dend.idxWholeCells(dend.indWholeCells);
wcAgglos = Agglo.fromSuperAgglo(dend.dendrites(dend.indWholeCells));

maxSegId = Seg.Global.getMaxSegId(param);
wcLUT = Agglo.buildLUT(maxSegId, wcAgglos);

somaWcIds = cellfun(@(ids) mode(nonzeros(wcLUT(ids))), somaAgglos);
assert(numel(somaWcIds) == numel(unique(somaWcIds)));

%% Completing and saving result
% TODO(amotta): Remove empty super-agglomerates

%% Debug inconsistencies
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
