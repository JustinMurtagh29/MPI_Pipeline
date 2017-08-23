%%
% script based on connectEM.axonQueryAnalysis by MB that collects
% axonsSmall agglomerates along the query flight paths
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

info = Util.runInfo(false);

nodeEvidenceStart = 13; % >
nodeEvidenceEnds = 53; % >

%% load data from MB

outputFolder = '/tmpscratch/mberning/axonQueryResults/';
Util.log('Loading query analysis result from MB from %s.', outputFolder);
resultMB = load([outputFolder 'postQueryState.mat']);
m = load('/gaba/scratch/mberning/axonQueryGeneration/axonsSmall.mat', ...
    'axonsSmall');
axonsSmall = m.axonsSmall;

%% rerun query overlap with axonsSmall

Util.log('Starting agglomeration overlap calculation.');
tic
[~, queryOverlap] = connectEM.queryAgglomerationOverlap( ...
    axonsSmall, resultMB.segmentsLeftover, resultMB.uniqueSegments, ...
    resultMB.neighboursStartNode);
t = toc;
Util.log('Finished agglomeration overlap in %f seconds', t);

eqClasses = [queryOverlap.ends.eqClasses];
occur = [queryOverlap.ends.occurences];
eqClasses = eqClasses(occur > nodeEvidenceEnds);
eqClasses = unique(eqClasses);

Util.log('%.3f of axonsSmall agglos were picked up by the queries.', ...
    length(eqClasses)./length(axonsSmall));

%% map for queries to axonsNew

Util.log('Getting query to axonsNew map.');

% only keep queries that have a start or ending
qIdx = ~resultMB.idxNoClearStart | ~resultMB.idxNoClearEnd;
queryOverlap.ends = queryOverlap.ends(qIdx);

% get the agglos that a query with either a clear start or end connected
startOrig = resultMB.queryOverlap.start(qIdx);
endOrig = resultMB.queryOverlap.ends(qIdx);
queryConn = arrayfun(@(x, y) ...
    cat(2, x.eqClasses(x.occurences > nodeEvidenceStart), ...
           y.eqClasses(y.occurences > nodeEvidenceEnds)), ...
    startOrig, endOrig, 'uni', 0);

% find to which axonsNew agglo a query belong by looking at the start axons
% agglo of the query and search for those ids in axonsNew
ids2AxonsNew = full(Seg.Global.eClassLookup(resultMB.eqClassCCfull));
tmp = cellfun(@(x)x(1), queryConn);
query2AxonsNew = ids2AxonsNew(tmp);

%% get the picked up agglomerates from the queries

Util.log('Determine agglos picked up on flight path.');

% get axonsSmall agglo idx and node evidence for corresponding queries
pickedUpAgglos = arrayfun(@(x)x.eqClasses(x.occurences > nodeEvidenceEnds)', ...
    queryOverlap.ends, 'uni', 0);
pickedUpEvidence = arrayfun(@(x)x.occurences(x.occurences > nodeEvidenceEnds)', ...
    queryOverlap.ends, 'uni', 0);

% assign agglos to queries with maximal node evidence (winner takes it all)
queryIdx = repelem((1:length(pickedUpAgglos))', ...
    cellfun(@length, pickedUpAgglos));
aggloIdx = cell2mat(pickedUpAgglos);
aggloEv = cell2mat(pickedUpEvidence);
s = sparse(queryIdx, aggloIdx, aggloEv);
[~, idx] = max(s, [], 1);
uAggloIdx = unique(aggloIdx);
aggloIdx2Query = idx(uAggloIdx)';
pickedUpAgglosWTA = accumarray(aggloIdx2Query, uAggloIdx, [], @(x){x});
pickedUpAgglosWTA = cellfun(@sort, pickedUpAgglosWTA, 'uni', 0);

% note that pickedUpAgglosWTA can be empty if all agglos were removed from
% one agglo (e.g. happens for the very last entry which is why the
% following line might be necessary)
pickedUpAgglosWTA(end+1:length(pickedUpAgglos)) = {[]};

%% add the axonsSmall agglos to the corresponding axonsNew agglos

Util.log('Updating axonsNew with picked up agglos.');

axonsNew = resultMB.axonsNew;
for i = 1:length(pickedUpAgglosWTA)
    if ~isempty(pickedUpAgglosWTA)
        axonsNew{query2AxonsNew(i)} = cat(1, ...
            axonsNew{query2AxonsNew(i)}, ...
            cell2mat(axonsSmall(pickedUpAgglosWTA{i})));
    end
end

% number of extended agglos
l_new = cellfun(@length, axonsNew);
l_old = cellfun(@length, resultMB.axonsNew);
numExtendedAgglos = sum(l_new > l_old);
% numExtendedAgglos should be the same as
% length(unique(query2AxonsNew(~cellfun(@isempty, pickedUpAgglosWTA))))

Util.log(['Extended %d agglos (%.3f of all agglos) with agglos ' ...
    'along flight path.'], ...
    numExtendedAgglos, numExtendedAgglos/length(axonsNew));

% number of picked up post query agglos
l = cellfun(@length, resultMB.axonsPostQuery);
tmp = cell2mat(resultMB.axonsPostQuery);
tmp2 = cell2mat(axonsNew);
pickedUpSegments = ismember(tmp, tmp2);
pickedUpSegments = mat2cell(pickedUpSegments, l, 1);
pickedUpAgglos = cellfun(@any, pickedUpSegments);

Util.log('Axons new picked up %.3f of axonsPostQuery agglos.', ...
    sum(pickedUpAgglos)./length(pickedUpAgglos));

%% save result

outFile = fullfile(outputFolder, 'axonQueryAnalysisBS.mat');
if ~exist(outFile, 'file')
    Util.log('Saving output to %s.', outFile);
    save(outFile', 'axonsNew', 'queryOverlap', 'info');
else
    warning(['The file %s already exists and will not be overwritten.' ...
        'Save the result manually.'], outFile);
end

Util.log('Finished picking up of agglomerates along flight path');

%% get boutons picked up by flight path
% this requires the result from the boutonAgglomerationScript

boutonFile = '/u/bstaffle/data/L4/Synapses/boutonAgglos.mat';
if ~exist(boutonFile, 'file')
    error('Run L4.Synapses.boutonAgglomerationScript first');
end
m = load(boutonFile, 'bAgglo');
bAgglo = m.bAgglo;
% combine overlapping agglos
bAgglo_ov = Seg.Global.combineEClasses(bAgglo);

Util.log('Starting bouton overlap calculation.');
tic
[partitionBoutons, queryOverlapBoutons] = connectEM.queryAgglomerationOverlap( ...
    bAgglo_ov, resultMB.segmentsLeftover, resultMB.uniqueSegments, ...
    resultMB.neighboursStartNode);
t = toc;
Util.log('Finished bouton overlap in %f seconds', t);

% save (append to file above if not yet in there)
m = load(outFile);
if isfield(m, 'queryOverlapBoutons')
    warning(['Bouton overlap already exists in file %s and will not be' ...
        ' overwritten.'], outFile)
else
    Util.log('Saving bouton overlap to %s.', outFile);
    m.queryOverlapBoutons = queryOverlapBoutons;
    m.partitionBoutons = partitionBoutons;
    m.bAgglo = bAgglo;
    save(outFile, '-struct', 'm');
end

% only keep queries that have a start or ending
qIdx = ~resultMB.idxNoClearStart | ~resultMB.idxNoClearEnd;
queryOverlapBoutons.ends = queryOverlapBoutons.ends(qIdx);

% get boutons from flight paths
pickedUpBoutons = arrayfun(@(x)x.eqClasses(x.occurences > nodeEvidenceEnds)', ...
    queryOverlapBoutons.ends, 'uni', 0);
flightPathBoutonIdx = cell2mat(pickedUpBoutons);

% get boutons in default post query axons
b2agglo_ov = L4.Synapses.getCollectedBoutons(resultMB.axonsNew, bAgglo_ov);

% all picked up boutons
pickedUpBoutonIdx = unique([flightPathBoutonIdx; find(b2agglo_ov > 0)]);

% statistics
stats.boutons.total = length(b2agglo_ov);
stats.boutons.picked_up_fraction = ...
    length(pickedUpBoutonIdx)/length(b2agglo_ov);
stats.boutons.picked_up_total = length(pickedUpBoutonIdx);

m = load(outFile);
if ~isfield(m, 'stats')
    m.stats = stats;
    save(outFile, '-struct', 'm');
end

%% write random agglos that changed to nml

l1 = cellfun(@length, axonsNew);
l2 = cellfun(@length, resultMB.axonsNew);

hasChanged = find(l1 > l2); % only stuff can be added

% rng('shuffle')
% idx = randi(length(hasChanged), 10, 1);
idx = [4808;1346;2933;7166;1643;5998;3412;6857;2509;6617];

m = load('/gaba/u/mberning/results/pipeline/20170217_ROI/allParameter.mat');
p = m.p;
m = load([p.saveFolder 'segmentMeta.mat'], 'point');
point = m.point';

agglos_new = axonsNew(hasChanged(idx));
agglos_old = resultMB.axonsNew(hasChanged(idx));
agglos = [agglos_new, agglos_old]';
agglos = agglos(:);

skel = L4.Agglo.agglo2Nml(agglos, point);

for i = 1:10
    skel.names{2*i - 1} = sprintf('AxonAgglo%d_flightPathExt', ...
        hasChanged(idx(i)));
    skel.names{2*i} = sprintf('AxonAgglo%d', hasChanged(idx(i)));
end

% skel.write('FlightPathSegmentCollection.nml')

%% wk mapping for agglos

tmpIdx = false(length(axonsSmall), 1);
tmpIdx(unique(aggloIdx)) = true;

m = load([p.saveFolder, 'globalSegSize.mat']);
maxSegId = length(m.segSize);

% create mapping classes:
% agglomerated >5um axons; axons small collected by flight path;
% uncollected axons small
components = cat(1, {vertcat(resultMB.axonsNew{:})}, ...
    {vertcat(axonsSmall{tmpIdx})}, ...
    {vertcat(axonsSmall{~tmpIdx})});
toZero = setdiff(0:maxSegId, cell2mat(components));
components = cat(1, toZero, components);
% WK.makeWKMapping( components, 'axonQueryFlightPathCollection' );

%% get volumes of collected/uncollected axons agglos

% total number of voxels
m = load([p.saveFolder 'segmentMeta.mat'], 'voxelCount');
voxelCount = m.voxelCount;
vxAxonsNew = sum(cellfun(@(x)sum(voxelCount(x)), axonsNew));
vxUncollected = sum(cellfun(@(x)sum(voxelCount(x)), axonsSmall(~tmpIdx)));
Util.log('Fraction of collected over all axon agglosvolume: %.3f', ...
    vxAxonsNew/(vxUncollected + vxAxonsNew));

% density wrt local cubes
axonsNewIds = cell2mat(axonsNew);
axonsNewCubeIds = Seg.Global.getCubeIds(p, axonsNewIds);
collectedSegPerCube = accumarray(axonsNewCubeIds, 1);
collectedVxPerCube = accumarray(axonsNewCubeIds, voxelCount(axonsNewIds));
collectedSegPerCube(end+1:numel(p.local)) = 0;
collectedVxPerCube(end+1:numel(p.local)) = 0;
uncollIds = cell2mat(axonsSmall(~tmpIdx));
uncollecteCubeIds = Seg.Global.getCubeIds(p, uncollIds);
uncollectedSegPerCube = accumarray(uncollecteCubeIds, 1);
uncollectedVxPerCube = accumarray(uncollecteCubeIds, voxelCount(uncollIds));
uncollectedSegPerCube(end+1:numel(p.local)) = 0;
uncollectedVxPerCube(end+1:numel(p.local)) = 0;

% reshape to a 2d representation
collectedSegPerCube = reshape(collectedSegPerCube, [11 17*13]);
uncollectedSegPerCube = reshape(uncollectedSegPerCube, [11 17*13]);
collectedVxPerCube = reshape(collectedVxPerCube, [11 17*13]);
uncollectedVxPerCube = reshape(uncollectedVxPerCube, [11 17*13]);

% % plot fraction of collected over all (visualized in 2d - z direction along
% % second axis)
% figure;
% a = subplot(2, 1, 1);
% imagesc(collectedSegPerCube./(collectedSegPerCube + uncollectedSegPerCube));
% title('Fraction of collected axon segments over all segments in axon agglo.')
% a.XTick = 1:17:221;
% colorbar
% Visualization.Figure.plotDefaultSettings()
% a = subplot(2, 1, 2);
% imagesc(collectedVxPerCube./(collectedVxPerCube + uncollectedVxPerCube));
% title('Fraction of collected axon voxels over all voxels in axon agglo.')
% a.XTick = 1:17:221;
% colorbar
% Visualization.Figure.plotDefaultSettings()