%% Settings
% Load parameter of pipeline run
load /tmpscratch/alik/pipelineRunApril2017/allParameter.mat;

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
%% Consolidate 3 redundant annotations into 1, correspondences needed
gtConsolidated = connectEM.consolidateRedundantAnnotations(p, gt);
%% Remove all nodes from skeletons that are not part of consensus (and regroup)
outputFolder=fullfile(newGTfolder,'consolidated');
if ~isdir(outputFolder)
    mkdir(outputFolder)
end
connectEM.restrictToConsolidationAndWriteAsSkeletons(gtConsolidated, p, pT, outputFolder);

%%Only run up until here 25.09.2017
%% First round of ground truth data proofreading (add leftover segments)
% Need to load "initial" instead of "afterCompletion" in first cell above
outputFolder = ['trainingDataOut' filesep];
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

%% Second round of ground truth data proofreading (inspect merger)
outputFolder = ['trainingDataOut2' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i=1:length(pT.local)
    for j=1:length(pT.local(i).trainFile)
        [~, originalFilename] = fileparts(pT.local(i).trainFile{j});
        skelFile = [outputFolder originalFilename '.nml'];
        mappingFile = [outputFolder originalFilename '.txt'];
        % Read skeleton in current state
        skel = skeleton(pT.local(i).trainFile{j});
        % Keep only trees in merged segments
        mergedSegments = gt(i,j).mergedSegments;
        treeIdx = ~cellfun(@(x)any(ismember(x, mergedSegments)), gt(i,j).segIdsOfGTnodes);
        skelNew = skel.deleteTrees(treeIdx);
        segIdsOfGTnodes = gt(i,j).segIdsOfGTnodes(~treeIdx);
        % Set node radius to 10 for all nodes
        for k=1:length(skelNew.nodesNumDataAll)
            skelNew.nodesNumDataAll{k}(:,2) = 10;
        end
        % Write comments to each node within each merged segment
        for k=1:length(mergedSegments)
            nodeIdx = cellfun(@(x)find(ismember(x, mergedSegments(k))), segIdsOfGTnodes, 'uni', 0);
            for n=1:length(nodeIdx)
                if ~isempty(nodeIdx{n})
                    % Mark one (random) node within each participating tree
                    % by setting comment and increasing node size
                    randIdx = randi(length(nodeIdx{n}), 1);
                    skelNew.nodesAsStruct{n}(nodeIdx{n}(randIdx)).comment = ...
                        ['merger' num2str(k, '%.3i')];
                    skelNew.nodesNumDataAll{n}(nodeIdx{n}(randIdx),2) = 30;
                end
            end
        end
        skelNew.write(skelFile);
        % Add mapping showing only merged segments
        script = WK.makeMappingScript(segMeta.maxSegId, num2cell(mergedSegments));
        fileHandle = fopen(mappingFile, 'w');
        fwrite(fileHandle, script);
        fclose(fileHandle);
    end
end

%% Third step: Read data after completion (execute first + second cell again)
% and then: modify according to merged segment decisions

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/mergedSegmentDecision/';
files = dir([newGTfolder '*.nml']);
mergerNumberOld = arrayfun(@(x)size(x.mergedSegments,1), gt);

% Determine annotated bounding box for each file (to match to training
% region below)
for i=1:length(files)
    skel = skeleton([newGTfolder files(i).name]);
    nodes = cat(1, skel.nodes{:});
    bbox{i} = [min(nodes(:,1:3),[],1); max(nodes(:,1:3),[],1)]';
    [~, pTbboxIdx(i)] = min(cellfun(@(x)norm((bbox{i} - x)'), {pT.local(:).bboxSmall}));
    uniqueComments = skel.getUniqueComments();
    temp = cellfun(@(x)regexp(x, 'merger(\d\d\d).+', 'tokens'), uniqueComments{:,1}, 'uni', 0);
    temp = temp(~cellfun(@isempty,temp));
    temp = cellfun(@(x)str2double(x{1}), temp);
    maxMergerId = max(temp(:));
    % Match to original tracing based on number of merger corrected
    % Next time find some better way to make sure relation is right, works here
    % because each file has different number of merger (TODO)
    [row, col] = find(mergerNumberOld == maxMergerId);
    pT.local(row).trainFileMerger{col} = [newGTfolder files(i).name];
end

% Cleanup, getting confused
clear row col files bbox nodes skel i j a ans temp taskIds tracer annotationIds newGTfolder maxMergerId pTbboxIdx mergerNumberOld uniqueComments;
save('/home/mberning/Desktop/tempMerger.mat');

%% Apply merger decision to segment lists and generate new skeletons
load('/home/mberning/Desktop/tempMerger.mat');
skelModified = connectEM.applyMergerDecision(pT);

% Then save for inspection again
outputFolder = ['trainingDataOut3' filesep];
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
for i=1:length(pT.local)
    for j=1:length(pT.local(i).trainFile)
        pT.local(i).trainFile{j} = [outputFolder 'trainingRegion' num2str(i) '_instance' num2str(j) '.nml'];
        skelModified(i,j) = skelModified(i,j).clearComments();
        skelModified(i,j) = skelModified(i,j).deleteEmptyTrees();
        skelModified(i,j).write(pT.local(i).trainFile{j});
    end
end
clear i j;

% Transitition to new pipeline run
gtNew = connectEM.getContinuityLabelsFromNml(p, pT);

% Consolidate 3 redundant annotations into 1
gtConsolidated = connectEM.consolidateRedundantAnnotations(p, gtNew);

% Remove all nodes from skeletons that are not part of consensus (and regroup)
connectEM.restrictToConsolidationAndWriteAsSkeletons(gtConsolidated, p, pT, outputFolder);

% See denseCubesExpertAnnotation for further steps taken on new pipeline
% run (more segments due to CC fix) based on nml written in last function above:
% https://gitlab.mpcdf.mpg.de/connectomics/pipeline/issues/32
