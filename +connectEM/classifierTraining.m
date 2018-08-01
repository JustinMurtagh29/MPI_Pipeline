%% Settings
rng default;
% Load parameter of pipeline run
load /gaba/u/mberning/results/pipeline/20170217_ROI/allParameterWithSynapses.mat;

% Extract information of dense cube annotations
newGTfolder = '+connectEM/trainingData/afterFNandMergerAnnotationMB3/';

% Define bounding box of training regions as passed to annotators
% +1 for wk vs. matlab coordinate system offset
pT.local(1).bboxSmall = [4133 3963 2253; 4578 4408 2432]'+1;
pT.local(2).bboxSmall = [4438 1320 893; 4883 1765 1072]'+1;
pT.local(3).bboxSmall = [1824 6673 1239; 2269 7118 1418]'+1;

pT.local(1).trainFile{1} = [newGTfolder 'region-1.nml'];
pT.local(2).trainFile{1} = [newGTfolder 'region-2.nml'];
pT.local(3).trainFile{1} = [newGTfolder 'region-3.nml'];

segMeta = load([p.saveFolder 'segmentMeta.mat']);
segMeta.point = segMeta.point';

%% Collect training data from segmentation + nmls
% Get a list of labels from each (redundant) training file in each region
gt = connectEM.getContinuityLabelsFromNml(p, pT);

%% Collect features from SynEM feature calculation for all edges in GT
% Note this will only keep borders > 10 voxel and sorts out edges that have
% multiple borders
gt = connectEM.getFeaturesForEdges(p, gt);
gt = rmfield(gt, {'segIdsOfGTnodes', 'segIdsGT', 'mergedSegments', 'leftSegments'});

%% Divide GT into training an test set (80%/20%) in each dense cube

train = struct();
test = struct();
for i=1:length(gt)
    idxTrain = randperm(numel(gt(i).labels), floor(numel(gt(i).labels).*0.8));
    idxTest = setdiff(1:numel(gt(i).labels), idxTrain);
    fieldNames = fieldnames(gt);
    for j=1:length(fieldNames)
        train(i).(fieldNames{j}) = gt(i).(fieldNames{j})(idxTrain,:);
        test(i).(fieldNames{j}) = gt(i).(fieldNames{j})(idxTest,:);
    end
end
train = Util.concatStructs(1, train(1), train(2), train(3));
test = Util.concatStructs(1, test(1), test(2), test(3));

% NOTE(amotta): At this point the label vector of the training set
% indicates not only positive and negative samples, but also contains zeros
% that make up ~15 % of the entires.
%   What's the meaning of these zeros? It's not clear to me what value
% these samples add for the training of a classifier. Shouldn't we just
% remove them from the training set?

%% Augment training & test set (by adding features for "inverted" edges)

load([p.saveFolder 'SynapseClassifier.bkp'], 'fm', '-mat');
fm.areaT = 10;
train.classFeaturesInv =  fm.invertDirection(train.classFeatures);
train.rawFeaturesInv =  fm.invertDirection(train.rawFeatures);
test.classFeaturesInv =  fm.invertDirection(test.classFeatures);
test.rawFeaturesInv =  fm.invertDirection(test.rawFeatures);

%% Train classifier
classifierNewFeatures = connectEM.trainClassifier( ...
    cat(1,cat(2, train.rawFeatures, train.classFeatures),cat(2, train.rawFeaturesInv, train.classFeaturesInv)), ...
    cat(1, train.labels, train.labels));

%% Determine whether compacted classifier works here
dateString = datestr(clock,30);
classifier = compact(classifierNewFeatures);
save(['/gaba/u/mberning/results/edgeClassifier/' dateString '.mat'], 'classifier');

a = classifier.predict(cat(2, test.rawFeatures, ...
    test.classFeatures));
b = classifierNewFeatures.predict(cat(2, test.rawFeatures, ...
    test.classFeatures));

all(a(:) == b(:))

%% Save everything for later checks
save(['/gaba/u/mberning/results/edgeClassifier/' dateString '_workspace.mat'], '-v7.3');

