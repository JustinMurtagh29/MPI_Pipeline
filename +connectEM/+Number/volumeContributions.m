% Calculates the contributions of nuclei, somata, and blood vessels to the
% total segmented volume.
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
connFile = fullfile(rootDir, 'connectomeState', 'connectome_axons-19-a-partiallySplit-v2_dendrites-wholeCells-03-v2-classified_SynapseAgglos-v8-classified.mat');

minSynCount = 10;
plotSynCountHists = true;

% See connectEM.WholeCell.Script.separateSomata
% 18c1bc8c43fb633af38dfa1ee60aba93424a9581
somaFile = fullfile(rootDir, 'aggloState', 'somata_07.mat');

vesselThreshold = 0.5;
nucleiThreshold = 0.5;

% EM data margin around processed bounding box
segEmFov = [25, 25, 10];

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

conn = connectEM.Connectome.load(param, connFile);

% Segment sizes
segSizes = Seg.Global.getSegToSizeMap(param);
voxelSize = param.raw.voxelSize / 1E3;
voxelVol = prod(voxelSize);

% Somata
somata = load(somaFile);
somata = somata.somata;
somata = Agglo.fromSuperAgglo(somata);
somata = unique(cell2mat(somata));

% Heurist nuclei and blood vessel detection
heuristics = fullfile(rootDir, 'heuristicResult.mat');
heuristics = load(heuristics);
heuristics.segIds = double(heuristics.segIds);

nuclei = heuristics.segIds(heuristics.nucleiScore > nucleiThreshold);
vessel = heuristics.segIds(heuristics.vesselScore > vesselThreshold);

neuropilSegIds = setdiff( ...
    1:numel(segSizes), cat(1, somata, nuclei, vessel));

axonSegIds = cell2mat(conn.axons);
dendAndSpineSegIds = cell2mat(conn.dendrites(~ismember( ...
    conn.denMeta.targetClass, {'AxonInitialSegment', 'Somata'})));

%% Dataset extent
rawSize = 1 + diff(param.bbox, 1, 2)';
rawSize = rawSize + 2 * segEmFov;
rawSizeUm = round(rawSize .* voxelSize, 1) %#ok

%% Calculate contributions
volT = struct;
volT(1).name = 'SegEM volume';
volT(1).vol = prod(diff(param.bbox, 1, 2) + 1);

volT(2).name = 'Watershed borders';
volT(2).vol = volT(1).vol - sum(segSizes);

volT(3).name = 'Somata';
volT(3).vol = sum(segSizes(somata));

volT(4).name = 'Nuclei';
volT(4).vol = sum(segSizes(nuclei));

volT(5).name = 'Blood vessels';
volT(5).vol = sum(segSizes(vessel));

volT(6).name = 'Neuropil';
volT(6).vol = sum(segSizes(neuropilSegIds));

volT = struct2table(volT);
volT.percent = 100 * volT.vol / volT.vol(1);
volT.vol = volT.vol * voxelVol;

% Calculate volume fractions
% without taking into account watershed borders
nonBorderVol = volT.vol(1) - volT.vol(2);
volT.borderCorrPercent = nan(size(volT.percent));
volT.borderCorrPercent(3:end) = 100 * volT.vol(3:end) ./ nonBorderVol;
disp(volT)

%% Volumetric axon to dendrite ratio
axonVol = sum(segSizes(axonSegIds));
dendAndSpineVol = sum(segSizes(dendAndSpineSegIds));
ourDendAndSpineToAxonVolumeRatio = dendAndSpineVol / axonVol %#ok

% From Braitenberg & Schüz. (1998) Cortex: Statistics and geometry of
% neuronal connectivity. Second thoroughly revised edition.
bsDendToAxonVolumeRatio = 35 / 34 %#ok
bsDendAndSpineToAxonVolumeRatio = (35 + 14) / 34 %#ok
