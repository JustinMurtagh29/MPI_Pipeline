% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');
asiFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-linearized_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified__20190227T082543_asiT.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_02_v3_auto.mat');

outDir = '/tmpscratch/amotta/l4/2019-10-09-axons-dendrites-connectome-in-hdf5';

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

voxelSize = param.raw.voxelSize;
maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);
segSizes = Seg.Global.getSegToSizeMap(param);

[conn, syn, axonClasses] = ...
    connectEM.Connectome.load(param, connFile);
[linConn, linSyn, linConnFile] = ...
    connectEM.Consistency.loadConnectome(param);

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

asiT = asiT(asiT.area > 0, :);
asiT = connectEM.Consistency.Calibration.apply(asiT);

%% Export axon and dendrite reconstructions
clear cur*;

% Axons
curOutFile = fullfile(outDir, 'axons.hdf5');
categoricalToHdf5(curOutFile, '/axons/class', conn.axonMeta.axonClass);
agglosToHdf5(curOutFile, '/axons/agglomerate', conn.axons);
superAgglosToHdf5(curOutFile, '/axons/skeleton', axons);
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

% Dendrites
curShVols = cellfun(@(ids) sum(segSizes(ids)), shAgglos);
curShVols = single(curShVols * prod(voxelSize / 1E3));

curOutFile = fullfile(outDir, 'dendrites_v2.hdf5');
arrayToHdf5(curOutFile, '/dendrites/neuronId', uint32(conn.denMeta.cellId));
categoricalToHdf5(curOutFile, '/dendrites/class', conn.denMeta.targetClass);
agglosToHdf5(curOutFile, '/dendrites/agglomerate', conn.dendrites);
superAgglosToHdf5(curOutFile, '/dendrites/skeleton', dendrites);
agglosToHdf5(curOutFile, '/spineHeads/agglomerate', shAgglos);
arrayToHdf5(curOutFile, '/spineHeads/volume', curShVols);
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

%% Export connectome
clear cur*;

curSynT = connectEM.Connectome.buildSynapseTable(conn, syn);
curLinSynT = connectEM.Connectome.buildSynapseTable(linConn, linSyn);
curShLUT = Agglo.buildLUT(maxSegId, shAgglos);

curCalcPos = @(type) uint32(round( ...
    connectEM.Synapse.calculatePositions(param, syn, type)));

curOut = table;
curOut.type = syn.synapses.type;
curOut.preSegIds = syn.synapses.presynId;
curOut.postSegIds = syn.synapses.postsynId;

curOut.prePosition = curCalcPos('preRecenter');
curOut.synPosition = curCalcPos('border');
curOut.postPosition = curCalcPos('postRecenter');

curOut.preAxonId(:) = uint32(0);
curOut.preAxonId(curSynT.id) = curSynT.preAggloId;

curOut.preSplitAxonId(:) = uint32(0);
curOut.preSplitAxonId(curLinSynT.id) = curLinSynT.preAggloId;

curOut.postDendriteId(:) = uint32(0);
curOut.postDendriteId(curSynT.id) = curSynT.postAggloId;

curOut.postSpineHeadId = uint32(cellfun( ...
    @(segIds) max(curShLUT(segIds)), ...
    curOut.postSegIds));

curOut.asiArea(:) = nan;
curOut.asiArea(asiT.id) = asiT.area;

curOutFile = fullfile(outDir, 'synapses_v2.hdf5');
curOut = table2struct(curOut, 'ToScalar', true);
structToHdf5(curOutFile, '/synapses', curOut);
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

%% Export segment positions
clear cur*;

curOutFile = fullfile(outDir, 'segments.hdf5');
numericToHdf5(curOutFile, '/segments/position', uint32(segPoints));
numericToHdf5(curOutFile, '/segments/voxelCount', uint32(segSizes));
infoToHdf5(curOutFile, info);
Util.protect(curOutFile);

%% Utilities
function superAgglosToHdf5(outFile, group, agglos)
    % Sanity checks
    assert(not(any(arrayfun(@(a) any(a.nodes(:, 4) <= 0), agglos))));
    assert(all(arrayfun(@(a) all(all(a.nodes(:, 1:3) >= 1)), agglos)));
    assert(all(arrayfun(@(a) all(all(a.edges > 0)), agglos)));
    
    agglos = rmfield(agglos, setdiff( ...
        fieldnames(agglos), {'nodes', 'edges'}));
    
    for curIdx = 1:numel(agglos)
        agglos(curIdx).segIds = reshape( ...
            uint32(agglos(curIdx).nodes(:, 4)), [], 1);
        agglos(curIdx).nodes = reshape( ...
            uint32(agglos(curIdx).nodes(:, 1:3)), [], 3);
        agglos(curIdx).edges = reshape(unique(sort( ...
            uint32(agglos(curIdx).edges), 2), 'rows'), [], 2);
        
        if isempty(agglos(curIdx).edges)
            % HACKHACKHACK(amotta): Due to some weird behaviour of MATLAB
            % or HDF5 it's not possible to store an empty edge list. Let's
            % insert a self-connection so that there's at least that...
            assert(not(isempty(agglos(curIdx).nodes)));
            agglos(curIdx).edges = ones(1, 2, 'uint32');
        end
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
