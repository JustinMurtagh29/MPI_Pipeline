%% Settings
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/afterHiwiConsolidation/';

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]'+1;
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]'+1;
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]'+1;

pT.local(1).trainFile{1} = [newGTfolder 'trainingRegion1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'trainingRegion2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'trainingRegion3.nml'];

segMeta = load([p.saveFolder 'segmentMeta.mat']);
segMeta.point = segMeta.point';
corr = Seg.Global.getGlobalCorrespondences(p);
clear newGTfolder;

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);

%% First round of ground truth data proofreading (add leftover segments)
% Need to load "initial" instead of "afterCompletion" in first cell above

outputFolder = ['trainingDataOut4' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i=1:length(pT.local)
    [~, originalFilename] = fileparts(pT.local(i).trainFile{1});
    mappingFileLeft = [outputFolder originalFilename 'left.txt'];
    mappingFileState = [outputFolder originalFilename 'state.txt'];
    % Add correspondences
    edges = cat(1, gt(i).edges, corr);
    % Generate skeleton from segIds in GT to avoid having multiple
    % nodes per segment and having edges as present in supervoxel graph
    connectEM.generateSkeletonFromAgglo(edges, segMeta.point, gt(i).segIdsGT, ...
        arrayfun(@(x)['Component' num2str(i) '_' num2str(x, '%.3i')], 1:length(gt(i).segIdsGT), 'uni', 0), ...
        outputFolder, segMeta.maxSegId);
    % Create one additional nml with left Todos (~300)
    skel = skeleton([outputFolder 'Component' num2str(i) '_001.nml'], 1);
    for k=1:length(gt(i).leftSegments)
        skel = skel.addTree(['Todo ' num2str(k, '%.4i')], segMeta.point(gt(i).leftSegments(k),:));
    end
    skel.write([outputFolder 'Component' num2str(i) '_left.nml'], 1:skel.numTrees, 1);
    % Add mapping showing only leftover segments
    script = WK.makeMappingScript(segMeta.maxSegId, num2cell(gt(i).leftSegments));
    fileHandle = fopen(mappingFileLeft, 'w+');
    fwrite(fileHandle, script);
    fclose(fileHandle);
    % And current state of agglomeration
    script = WK.makeMappingScript(segMeta.maxSegId, gt(i).segIdsGT);
    fileHandle = fopen(mappingFileState, 'w+');
    fwrite(fileHandle, script);
    fclose(fileHandle);
end

%% Settings: 2nd round for checking results
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/afterExpertAnnotation1/';

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]'+1;
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]'+1;
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]'+1;

pT.local(1).trainFile{1} = [newGTfolder 'trainingRegion1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'trainingRegion2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'trainingRegion3.nml'];

segMeta = load([p.saveFolder 'segmentMeta.mat']);
segMeta.point = segMeta.point';
corr = Seg.Global.getGlobalCorrespondences(p);
clear newGTfolder;

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT, 0);

%% First round of ground truth data proofreading (add leftover segments)
% Need to load "initial" instead of "afterCompletion" in first cell above

outputFolder = ['trainingDataOut5' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i=1:length(pT.local)
    [~, originalFilename] = fileparts(pT.local(i).trainFile{1});
    mappingFileLeft = [outputFolder originalFilename 'left.txt'];
    mappingFileState = [outputFolder originalFilename 'state.txt'];
    % Add mapping showing only leftover segments
    script = WK.makeMappingScript(segMeta.maxSegId, num2cell(gt(i).leftSegments));
    fileHandle = fopen(mappingFileLeft, 'w+');
    fwrite(fileHandle, script);
    fclose(fileHandle);
    % And current state of agglomeration
    script = WK.makeMappingScript(segMeta.maxSegId, gt(i).segIdsGT);
    fileHandle = fopen(mappingFileState, 'w+');
    fwrite(fileHandle, script);
    fclose(fileHandle);
end

%% Settings: 3rd round for checking results
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/afterExpertAnnotation6/';

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]'+1;
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]'+1;
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]'+1;

pT.local(1).trainFile{1} = [newGTfolder 'trainingRegion1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'trainingRegion2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'trainingRegion3.nml'];

segMeta = load([p.saveFolder 'segmentMeta.mat']);
segMeta.point = segMeta.point';
corr = Seg.Global.getGlobalCorrespondences(p);
clear newGTfolder;

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);

%% First round of ground truth data proofreading (add leftover segments)
% Need to load "initial" instead of "afterCompletion" in first cell above

outputFolder = ['trainingDataOut9' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i=1:length(pT.local)
    [~, originalFilename] = fileparts(pT.local(i).trainFile{1});
    skelFile = [outputFolder originalFilename '.nml'];
    mappingFileLeft = [outputFolder originalFilename 'left.txt'];
    mappingFileState = [outputFolder originalFilename 'state.txt'];
    % Generate skeleton from segIds in GT
    nodes = cellfun(@(x)segMeta.point(x,1:3), gt(i).segIdsGT, 'uni', 0);
    treeNames = arrayfun(@(x)['Component' num2str(i) '_' num2str(x, '%.3i')], 1:length(gt(i).segIdsGT), 'uni', 0);
    connectEM.generateSkeletonFromNodes(skelFile, nodes, treeNames);
    % Create one additional nml with left Todos (~300)
    if ~isempty(gt(i).leftSegments)
        skel = skeleton(pT.local(i).trainFile{1}, 1);
        for k=1:length(gt(i).leftSegments)
            skel = skel.addTree(['Todo ' num2str(k, '%.4i')], segMeta.point(gt(i).leftSegments(k),:));
        end
        skel.write([outputFolder 'Component' num2str(i) '_left.nml'], 1:skel.numTrees, 1);
    end
    % Add mapping showing only leftover segments
    script = WK.makeMappingScript(segMeta.maxSegId, num2cell(gt(i).leftSegments));
    fileHandle = fopen(mappingFileLeft, 'w+');
    fwrite(fileHandle, script);
    fclose(fileHandle);
    % And current state of agglomeration
    script = WK.makeMappingScript(segMeta.maxSegId, gt(i).segIdsGT);
    fileHandle = fopen(mappingFileState, 'w+');
    fwrite(fileHandle, script);
    fclose(fileHandle);
end

%% Settings: 4th round for checking results
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/afterTreeAnnotationByAM/';

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]'+1;
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]'+1;
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]'+1;

pT.local(1).trainFile{1} = [newGTfolder 'region-1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'region-2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'region-3.nml'];

segMeta = load([p.saveFolder 'segmentMeta.mat'], 'point', 'maxSegId');
segMeta.point = segMeta.point';

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);

%% Find edges where labels and edges disagree & write to skeleton file
outputFolder = ['+connectEM' filesep 'trainingData' filesep 'trainingDataOut11' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for i=1:length(gt)
    idx = find(gt(i).labels == -1 & gt(i).prob > .9);
    cc = mat2cell(gt(i).edges(idx,:), ones(size(gt(i).edges(idx,:),1),1), 2);
    treeNames = arrayfun(@(x)['region' num2str(i) '_fnCanidates_score' num2str(gt(i).prob(x), '%3.2f') '_component' num2str(x, '%.2i') ], idx, 'uni', 0);
    connectEM.generateSkeletonFromAgglo(gt(i).edges(idx,:), segMeta.point, cc, treeNames, ... 
        outputFolder, segMeta.maxSegId);

    idx = find(gt(i).labels == 1 & gt(i).prob < .1);
    cc = mat2cell(gt(i).edges(idx,:), ones(size(gt(i).edges(idx,:),1),1), 2);
    treeNames = arrayfun(@(x)['region' num2str(i) '_fpCanidates_score' num2str(gt(i).prob(x), '%3.2f') '_component' num2str(x, '%.2i') ], idx, 'uni', 0);
    connectEM.generateSkeletonFromAgglo(gt(i).edges(idx,:), segMeta.point, cc, treeNames, ... 
        outputFolder, segMeta.maxSegId);
end

