% Look at random examples of endings which were created by the
% splitting of a chiasma. Questions to answer:
% 
%   1) How many of these are there?
%   2) How many of these are splits?
%   3) How can we query them?
%      a) Seeding: Position and direction?
%      b) Attachment
%
% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;
rootDir = '/gaba/u/mberning/results/pipeline/20170217_ROI';
outputDir = '/home/amotta/Desktop/chiasma-healing';

%% loading parameters
paramFile = fullfile(rootDir, 'allParameter.mat');
param = load(paramFile, 'p');
param = param.p;

%% loading axons
axonFile = fullfile(rootDir, 'aggloState', 'axons_13_a.mat');
axons = load(axonFile, 'axons', 'indBigAxons');

axonIds = find(axons.indBigAxons);
axons = axons.axons(axonIds);

%% collecting and counting endings
endings = connectEM.Chiasma.Heal.buildTaskDefs(param, axons);

numEndings = size(endings, 1);
fprintf('Total number of endings: %d\n', numEndings);

%% exporting random examples
numExamples = 20;

rng(0);
randIds = randperm(size(endings, 1), numExamples);
randEndings = endings(randIds, :);

mkdir(outputDir);

for curIdx = 1:numExamples
    curSkel = skeleton();
    
    curAxonId = randEndings.axonId(curIdx);
    curNodeId = randEndings.nodeId(curIdx);
    
    curAxon = axons(curAxonId);
    curNodeCount = size(curAxon.nodes, 1);
    
    curComments = repmat({''}, [curNodeCount, 1]);
    curComments{curNodeId} = 'Ending';
    
    curSkel = curSkel.addTree( ...
        'Axon', curAxon.nodes(:, 1:3), ...
        curAxon.edges, [], [], curComments);
    curSkel = curSkel.addBranchpoint(curNodeId);
    curSkel = Skeleton.setParams4Pipeline(curSkel, param);
    
    curQueryStart = randEndings.position(curIdx, :);
    curQueryEnd = randEndings.direction(curIdx, :);
    curQueryEnd = curQueryEnd ./ norm(curQueryEnd);
    curQueryEnd = curQueryStart + 10 .* curQueryEnd;
    
    curSkel = curSkel.addTree( ...
        'Query', [curQueryStart; curQueryEnd], [1, 2], [0, 0, 1, 1]);
    
    curSkelFile = sprintf( ...
        '%d__axon-%d__node-%d.nml', ....
        curIdx, curAxonId, curNodeId);
    curSkel.write(fullfile(outputDir, curSkelFile));
end