%% Settings
% Load parameter of pipeline run
load /tmpscratch/alik/pipelineRunApril2017/allParameter.mat;
mainFolder='/tmpscratch/alik/pipelineRunApril2017/trainingData/';
% Extract information of dense cube annotations
newGTfolder = '/tmpscratch/alik/pipelineRunApril2017/trainingData/SecondRoundHiwiTracing/';
files = dir([newGTfolder '*.nml']);

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = Util.convertWebknossosToMatlabBbox([3781, 1525, 3214, 417, 417, 167]);%contains myelin
pT.local(2).bboxSmall = Util.convertWebknossosToMatlabBbox([215, 6630, 3736, 417, 417, 167]);
pT.local(3).bboxSmall = Util.convertWebknossosToMatlabBbox([4248, 249, 1050, 417, 417, 167]);%contains soma

% Determine annotated bounding box for each file (to match to training
% region below)
for i=1:9
    skel = skeleton([newGTfolder files(i).name]);
    nodes = cat(1, skel.nodes{:});
    bbox{i} = [min(nodes(:,1:3),[],1); max(nodes(:,1:3),[],1)]';
    [~, pTbboxIdx(i)] = min(cellfun(@(x)norm((bbox{i} - x)'), {pT.local(:).bboxSmall}));
end

files = {files(:).name};

pT.local(1).trainFile = cellfun(@(x)[newGTfolder x], files(pTbboxIdx == 1), 'uni', 0);
pT.local(2).trainFile = cellfun(@(x)[newGTfolder x], files(pTbboxIdx == 2), 'uni', 0);
pT.local(3).trainFile = cellfun(@(x)[newGTfolder x], files(pTbboxIdx == 3), 'uni', 0);

segMeta = load([p.saveFolder 'segmentMeta.mat']);
segMeta.point = segMeta.point';

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);
clear files bbox nodes skel i newGTfolder maxMergerId pTbboxIdx;

%% First round of ground truth data proofreading (add leftover segments)
% Need to load "initial" instead of "afterCompletion" in first cell above
outputFolder = [fullfile(mainFolder,'trainingDataOut') filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i=1:length(pT.local)
    for j=1:length(pT.local(i).trainFile)
        [~, originalFilename] = fileparts(pT.local(i).trainFile{j});
        skelFile = [outputFolder originalFilename '.nml'];
        mappingFile = [outputFolder originalFilename '.txt'];
        % Generate skeleton with additional nodes for not yet traced
        % segments
        skel = skeleton(pT.local(i).trainFile{j});
        % Set all old (already annoatated) skeletons to blue
        skel.colors = cellfun(@(x)[0 0 1 1], skel.colors, 'uni', 0);
        leftSegments = gt(i,j).leftSegments;
        for k=1:length(leftSegments)
            skel = skel.addTree(['Todo ' num2str(k, '%.4i')], segMeta.point(leftSegments(k),:));
        end
        skel.write(skelFile);
        % Add mapping showing only leftover segments
        script = WK.makeMappingScript(segMeta.maxSegId, num2cell(leftSegments));
        fileHandle = fopen(mappingFile, 'w');
        fwrite(fileHandle, script);
        fclose(fileHandle);
    end
end
%%Up to here so far for the first round of tracings

%% Consolidate 3 redundant annotations into 1, correspondences needed
gtConsolidated = connectEM.consolidateRedundantAnnotations(p, gt);
%% Remove all nodes from skeletons that are not part of consensus (and regroup)
outputFolder=fullfile(newGTfolder,'consolidated');
if ~isdir(outputFolder)
    mkdir(outputFolder)
end
connectEM.restrictToConsolidationAndWriteAsSkeletons(gtConsolidated, p, pT, outputFolder);

%% From here on comes from the expert annotation
newGTfolder = '/tmpscratch/alik/pipelineRunApril2017/trainingData/afterCompletion/consolidated/';
hiwiConsFolder=[fullfile('/tmpscratch/alik/pipelineRunApril2017/trainingData/','Round1ForHiwiConsolidation') filesep];
%Create the output folder
if ~exist(hiwiConsFolder, 'dir')
    mkdir(hiwiConsFolder);
end
%Load the correspondences
corr = Seg.Global.getGlobalCorrespondences(p);
pT.local=rmfield(pT.local,'trainFile');
pT.local(1).trainFile{1} = [newGTfolder 'trainingRegion1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'trainingRegion2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'trainingRegion3.nml'];
%% Collect training data from segmentation + nmls for the first round of giving back to Hiwis
% Get a list of labels from each (redundant) training file in each region
clear gt
gt = connectEM.getContinuityLabelsFromNml(p, pT);

for i=1:length(pT.local)
    [~, originalFilename] = fileparts(pT.local(i).trainFile{1});
    mappingFileLeft = [hiwiConsFolder originalFilename 'left.txt'];
    mappingFileState = [hiwiConsFolder originalFilename 'state.txt'];
    % Add correspondences
    edges = cat(1, gt(i).edges, corr);
    % Generate skeleton from segIds in GT to avoid having multiple
    % nodes per segment and having edges as present in supervoxel graph
    connectEM.generateSkeletonFromAgglo(edges, segMeta.point, gt(i).segIdsGT, ...
        arrayfun(@(x)['Component' num2str(i) '_' num2str(x, '%.3i')], 1:length(gt(i).segIdsGT), 'uni', 0), ...
        hiwiConsFolder, segMeta.maxSegId);
    % Create one additional nml with left Todos (~300)
    skel = skeleton([hiwiConsFolder 'Component' num2str(i) '_001.nml'], 1);
    for k=1:length(gt(i).leftSegments)
        skel = skel.addTree(['Todo ' num2str(k, '%.4i')], segMeta.point(gt(i).leftSegments(k),:));
    end
    skel.write([hiwiConsFolder 'Component' num2str(i) '_left.nml'], 1:skel.numTrees, 1);
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
