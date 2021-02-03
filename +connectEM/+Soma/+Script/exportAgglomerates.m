% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = fullfile(rootDir, 'allParameter.mat');
param = Util.load(param, 'p');

[conn, syn, axonClasses] = connectEM.Connectome.load(param, connFile);
dendrites = Util.load(conn.info.param.dendriteFile, 'dendrites');

%% Statistics
pdMask = conn.denMeta.targetClass == 'WholeCell';
pdInConnMask = pdMask & (conn.denMeta.synCount >= 10);
fprintf('Number of proximal dendrites: %d\n', sum(pdMask));
fprintf('Number of proximal dendrites with >= 10 synapses: %d\n', sum(pdInConnMask));

somaMask = conn.denMeta.targetClass == 'Somata';
somaHasPdInConn = conn.denMeta.cellId(pdInConnMask);
somaHasPdInConn = ismember(conn.denMeta.cellId(somaMask), somaHasPdInConn);
fprintf('Number of soma agglomerates: %d\n', sum(somaMask));

%% Export all somata
clear cur*;
curVols = Seg.Global.getSegToSizeMap(param);

somaAgglos = conn.denMeta.parentId(somaMask);
somaAgglos = dendrites(somaAgglos);

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

curDigits = ceil(log10(1 + numel(somaAgglos)));
for curIdx = 1:numel(somaAgglos)
    curSkel = somaAgglos(curIdx);
    curSkel = SuperAgglo.connect(curSkel);
    
    curMask = not(isnan(curSkel.nodes(:, 4)));
    curPos = curVols(curSkel.nodes(curMask, 4));
    curPos = sum((curPos ./ sum(curPos)) .* curSkel.nodes(curMask, 1:3));
    
    curYesNo = {'no', 'yes'};
    curHasPdInConn = somaHasPdInConn(curIdx);
    curHasPdInConn = curYesNo{1 + curHasPdInConn};
    
    curName = sprintf( ...
        'Soma %0*d. PDs in connectome: %s. Center of mass: %d, %d, %d', ...
        curDigits, curIdx, curHasPdInConn, round(curPos));
    skel = skel.addTree(curName, curSkel.nodes(:, 1:3), curSkel.edges);
end
