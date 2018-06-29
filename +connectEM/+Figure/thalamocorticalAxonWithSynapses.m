% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a_dendrites-wholeCells-03-v2-classified_spine-syn-clust.mat');
shFile = fullfile(rootDir, 'aggloState', 'dendrites_wholeCells_01_spine_attachment.mat');

axonIds = 7680;
outputDir = '/home/amotta/Desktop';

info = Util.runInfo();

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

maxSegId = Seg.Global.getMaxSegId(param);
segPoints = Seg.Global.getSegToPointMap(param);

[conn, syn] = connectEM.Connectome.load(param, connFile);

shAgglos = load(shFile, 'shAgglos');
shAgglos = shAgglos.shAgglos;

axons = load(conn.info.param.axonFile);
axons = axons.axons;

%% Prepare synapse list
shLUT = Agglo.buildLUT(maxSegId, shAgglos);
synT = connectEM.Connectome.buildSynapseTable(conn, syn);

synT.shId = cellfun( ...
    @(segIds) max(shLUT(segIds)), ...
    syn.synapses.postsynId(synT.id));
synT(~synT.shId, :) = [];

%% Export to webKNOSSOS
numDigits = ceil(log10(1 + numel(axonIds)));
uncell = @(c) c{1};

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = skel.setDescription(sprintf( ...
    '%s (%s)', info.filename, info.git_repos{1}.hash));

for curIdx = 1:numel(axonIds)
    curId = axonIds(curIdx);
    
    curAxon = axons(conn.axonMeta.parentId(curId));
   [~, curFlights] = Superagglos.splitIntoAgglosAndFlights(curAxon);
   
    curFlights = uncell(curFlights);
    curFlightCount = numel(curFlights);
    
    curAxon = rmfield(curAxon, setdiff( ...
        fieldnames(curAxon), {'nodes', 'edges'}));
    curSkel = Superagglos.toSkel(cat(1, curAxon, curFlights), skel);
    
    curSkel.names{1} = sprintf('Axon %d', curId);
    curSkel.colors{1} = [0, 0, 1, 1];
    
    curSkel.names((end - curFlightCount + 1):end) = arrayfun( ...
        @(flightId) sprintf('Flight %d', flightId), ...
        1:curFlightCount, 'UniformOutput', false);
    curSkel.colors((end - curFlightCount + 1):end) = {[1, 0, 0, 1]};
    
    curSynT = synT(synT.preAggloId == curId, :);
    curSyns = shAgglos(curSynT.shId);
    
    curNodes = cellfun( ...
        @(segIds) segPoints(segIds, :), ...
        curSyns, 'UniformOutput', false);
    curSkel = Skeleton.fromMST( ...
        curNodes, param.raw.voxelSize, curSkel);
    
    curSkel.names((end - height(curSynT) + 1):end) = arrayfun( ...
        @(id) sprintf('Spine head %d', id), ...
        curSynT.shId, 'UniformOutput', false);
    curSkel.colors((end - height(curSynT) + 1):end) = {[1, 1, 0, 1]};
    
    curSkelName = sprintf( ...
        '%0*d_axon-%d.nml', numDigits, curIdx, curId);
    curSkel.write(fullfile(outputDir, curSkelName));
end
