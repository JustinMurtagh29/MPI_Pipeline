% Written by
%   Alessandro Motta <alessandro.motta@brain.mpg.de>
clear;

%% Configuration
rootDir = '/tmpscratch/sahilloo/data/Mk1_F6_JS_SubI_v1/pipelineRun_mr2e_wsmrnet/';
param = load(fullfile(rootDir, 'allParameter.mat'));
param = param.p;

nmlDir = fullfile(param.saveFolder, 'tracings', 'spine-attachment');
dendFile = fullfile(param.saveFolder, 'aggloMat', '20191227T134319_ga_20191224T235355optimParams_agglomeration/20191227T220548_results_auto-spines_v1.mat');
skelFile = fullfile(param.saveFolder, 'spineAttachment', 'debug-spine-attachment.nml');

exportCount = 200;

info = Util.runInfo();
Util.showRunInfo(info);

%% Loading data
points = Seg.Global.getSegToPointMap(param);
dend = load(dendFile, 'shAgglos', 'shEdges', 'edges', 'attached');

%% Export random subset of spines
clear cur*;

rng(0);
randIds = numel(dend.shAgglos);
randIds = randperm(randIds, exportCount);

numDigits = ceil(log10(1 + numel(randIds)));

skel = skeleton();
skel = Skeleton.setParams4Pipeline(skel, param);
skel = Skeleton.setDescriptionFromRunInfo(skel, info);

for curIdx = 1:numel(randIds)
    curId = randIds(curIdx);
    curAttached = dend.attached(curId);
    
    curName = sprintf( ...
        '%0*d. Spine head %d. Attached to %d.', ...
        numDigits, curIdx, curId, curAttached);
    
    curShSegIds = double(dend.shAgglos{curId});
    
    % HACKHACKHACK(amotta): Prior to
    % 
    % commit ce48b628219cd5f2924bc96ded6efa1e669fabfe
    % Author: Alessandro Motta <alessandro.motta@brain.mpg.de>
    % Date:   Thu Dec 12 13:41:56 2019 +0100
    % 
    %     Save both edges in spine head agglomerate and along neck
    %
    % both the `edges` and `shEdges` contained the neck edges. The edges
    % within the spine head were lost. So, let's generate some random edges
    % to ensure connectivity for now.
    curShEdges = horzcat( ...
        reshape(curShSegIds(1:(end - 1)), [], 1), ...
        reshape(curShSegIds(2:end), [], 1));
    
    curEdges = dend.shEdges{curId};
    curEdges = curEdges(all(curEdges, 2), :);
    curEdges = vertcat(curShEdges, curEdges); %#ok
    
    curSegIds = union(curShSegIds, curEdges);
    curNodes = points(curSegIds, :);
    
   [~, curEdges] = ismember(curEdges, curSegIds);
    assert(all(curEdges(:)));
    
    % remove duplicate edges
    curEdges = unique(curEdges, 'rows');
    skel = skel.addTree(curName, curNodes, curEdges);
end

skel.write(skelFile);
