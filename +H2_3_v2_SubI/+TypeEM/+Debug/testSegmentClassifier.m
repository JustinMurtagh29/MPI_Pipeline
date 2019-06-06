% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
config = 'ex144_08x2_mrNet';
classifierFile = '/tmpscratch/amotta/l4/2018-10-10-mrnet-pipeline-run/segmentClassifier.mat';

gtNmlFile = Util.getGitReposOnPath();
gtNmlFile = fullfile( ...
    gtNmlFile{1}, 'data', 'tracings', ...
    'ex144-08x2-mrNet', 'spine-heads-1.nml');

info = Util.runInfo();
Util.showRunInfo(info);

%% Load data
config = loadConfig(config);
param = config.param;

classifiers = load(classifierFile);

%% Load ground truth
clear cur*;

curNml = slurpNml(gtNmlFile);
curNodes = NML.buildNodeTable(curNml);

curNodes.coord = curNodes.coord + 1;
curNodes.segId = Seg.Global.getSegIds(param, curNodes.coord);
assert(all(curNodes.segId));

curBox = curNml.parameters.userBoundingBox;
curBox = { ...
    curBox.topLeftX, curBox.topLeftY, curBox.topLeftZ, ...
    curBox.width, curBox.height, curBox.depth};
curBox = cellfun(@str2double, curBox);

curBox = Util.convertWebknossosToMatlabBbox(curBox);
curseg = loadSegDataGlobal(param.seg, curBox);

posSegIds = reshape(unique(curNodes.segId), [], 1);
negSegIds = reshape(setdiff(curseg, [0; posSegIds]), [], 1);

%% Build ground truth structure
gt = struct;
gt.class = classifiers.classes;
gt.segId = cat(1, double(posSegIds), double(negSegIds));
gt.label = zeros(numel(gt.segId), numel(gt.class));

curMask = gt.class == 'spinehead';
gt.label(1:numel(posSegIds), curMask) = +1;
gt.label((numel(posSegIds) + 1):end, curMask) = -1;

%% Load features
gt = TypeEM.GroundTruth.loadFeatures( ...
    param, classifiers.featureSetName, gt);
gt.scores = TypeEM.Classifier.apply(classifiers, gt.featMat);

%% Plot precision / recall
clear cur*;
curLimits = [0, 100];

[~, fig] = TypeEM.Classifier.evaluate(param, classifiers, gt);

fig.Color = 'white';
fig.Position(3:4) = [410, 390];

ax = findobj(fig.Children, 'Type', 'Axes');
leg = findobj(fig.Children, 'Type', 'Legend');

axis(ax, 'square');
ax.TickDir = 'out';
xlim(ax, curLimits);
ylim(ax, curLimits);

leg.Box = 'off';
leg.Location = 'SouthWest';

title(ax, ...
    {info.filename; info.git_repos{1}.hash}, ...
    'FontWeight', 'normal', 'FontSize', 10);

%% Look at false positives
clear cur*;

curCount = 50;
curDigits = ceil(log10(1 + curCount));
curPoints = Seg.Global.getSegToPointMap(param);

curMask = gt.class == 'spinehead';

fp = table;
fp.segId = gt.segId;
fp.label = gt.label(:, curMask);
fp.score = gt.scores(:, curMask);

fp = fp(fp.label < 0, :);
fp = sortrows(fp, 'score', 'descend');

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

curSegIds = fp.segId(1:curCount);
curScores = fp.score(1:curCount);

curNodes = curPoints(curSegIds, :);
curNames = arrayfun( ...
    @(idx, segId, score) sprintf( ...
        '%0*d. Segment %d. Score %.3f', ...
        curDigits, idx, segId, score), ...
	(1:curCount)', curSegIds, curScores, ...
    'UniformOutput', false);

skel = skel.addNodesAsTrees(curNodes, curNames);
% skel.write('/home/amotta/Desktop/spine-head-fps.nml');
