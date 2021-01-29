% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outDir = '/tmpscratch/amotta/l4/2021-01-22-export-of-flight-paths';
axonFile = fullfile(rootDir, 'aggloState', 'axons_19_a_partiallySplit_v2.mat');

minNodeDistNm = 250;
runId = datestr(now, 30);

info = Util.runInfo();
Util.showRunInfo(info);

%% Initialization
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

outDir = fullfile(outDir, runId);
assert(not(exist(outDir, 'dir')));
assert(mkdir(outDir));

%% Extract flight paths of axons
axons = Util.load(axonFile, 'axons');

Util.log('Extracting flight paths');
[~, axonFlights] = Superagglos.splitIntoAgglosAndFlights(axons);
Util.log('Done');

axonIds = cellfun(@numel, axonFlights);
axonIds = repelem(1:numel(axonFlights), axonIds);
axonIds = reshape(axonIds, [], 1);

axonFlights = cat(2, axonFlights{:});
axonFlights = axonFlights(:);

out = struct;
out.info = info;
out.axonIds = axonIds;
out.axonFlights = axonFlights;

outFile = fullfile(outDir, 'flightpaths.mat');
Util.saveStruct(outFile, out);
Util.protect(outFile);
Util.clear(out);

%% (Try to) build NML file with skeletons
axonDigits = ceil(log10(1 + numel(axons)));
flightDigits = accumarray(axonIds(:), 1, size(axons(:)));
flightDigits = ceil(log10(max(flightDigits) + 1));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

axonId = -1;
flightId = -1;
for curIdx = 1:numel(axonFlights)
    curAxonId = axonIds(curIdx);
    
    if axonId ~= curAxonId
        axonId = curAxonId;
        flightId = 0;
    end
    
    curFlight = axonFlights(curIdx);
    % NOTE(amotta): Skip empty flight paths
    if isempty(curFlight.edges); continue; end
    flightId = flightId + 1;
    
    curFlight = Skeleton.reduceNodeDensity( ...
        curFlight, minNodeDistNm, param.raw.voxelSize);
    
    curName = sprintf( ...
        'Axon %0*d, flight %0*d', ...
        axonDigits, axonId, ...
        flightDigits, flightId);
    skel = skel.addTree(curName, ...
        curFlight.nodes(:, 1:3), ...
        curFlight.edges);
end

outFile = fullfile(outDir, 'flightpaths-axons.nml');
skel.write(outFile);
Util.protect(outFile);

