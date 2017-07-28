function trainingPipeline( dataFolder, saveFolder, fm, noTest, ...
    labelType, ensembleArgs, cluster )
%TRAININGPIPELINE Training pipeline for SynEM classifiers.
% INPUT dataFolder: string
%           Path to training data folder.
%       saveFolder: string
%           Path to folder where features are saved.
%       featureMap: Interface.FeatureMap object
%       noTest: (Optional) logical
%           Do not evaluate the classifier on the test set.
%           (Default: false)
%       labelType: (Optional) string
%           See 'type' in SynEM.Util.getTrainingDataFrom.
%           (Default: 'direction')
%       ensembleArgs: (Optional) [2Nx1] cell
%           Cell array of name value pairs for
%           SynEM.Classifier.BoostedEnsemble.train.
%       cluster: (Optional) parallel.cluster object
%           Cluster object used for calculation.
%           (Default: getCluster('cpu') if this function exists)
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

if ~exist('noTest','var') || isempty(noTest)
    noTest = false;
end
if ~exist('labelType','var') || isempty(labelType)
    labelType = 'direction';
end
if ~exist('cluster','var') || isempty(cluster)
    cluster = [];
end

dataFolder = SynEM.Util.addFilesep(dataFolder);
saveFolder = SynEM.Util.addFilesep(saveFolder);

if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end
fprintf('[%s] Results are stored in %s.\n', datestr(now), saveFolder);

fprintf('[%s] Calculating features.\n', datestr(now));
job = SynEM.Training.calculateFeaturesForTrainingData( dataFolder, ...
    saveFolder, fm, true, cluster);

wait(job);

SynEM.Util.groupTrainingData(saveFolder);

fprintf('[%s] Training classifier.\n', datestr(now));
[X, y] = SynEM.Util.getTrainingDataFrom([saveFolder, 'train'], labelType);

if exist('ensembleArgs','var') && ~isempty(ensembleArgs)
    classifier = SynEM.Classifier.BoostedEnsemble.train(X, y, ...
        ensembleArgs{:});
else
    classifier = SynEM.Classifier.BoostedEnsemble.train(X, y);
end
classifier.options.fm = fm;

fprintf('[%s] Saving classifier.\n', datestr(now));
mkdir([saveFolder, 'classifier']);
save([saveFolder, 'classifier/LogitBoostFull.mat'],'classifier','-v7.3');
classifier = compact(classifier);
classifier = classifier.calculatePredVar();
save([saveFolder, 'classifier/Classifier.mat'],'classifier');

fprintf('[%s] Validation performance.\n', datestr(now));
if strcmp(labelType, 'augment') %use single labels for validation
    labelType = 'single';
end
[X, ~] = SynEM.Util.getTrainingDataFrom([saveFolder, 'val'], labelType);
[~,scoresVal] = classifier.predict(X); %#ok<ASGLU>

if noTest
    save([saveFolder, 'classifier/result.mat'], 'scoresVal');
else
    m = load('allParameter20141007.mat');
    p = m.p;
    fprintf('[%s] Test performance.\n', datestr(now));
    scoresTest = SynEM.Seg.predictCube(p, 67, fm, classifier); %#ok<NASGU>
    save([saveFolder, 'classifier/result.mat'], 'scoresTest', 'scoresVal');
end

%calculate prediction classifier
imp = classifier.ens.predictorImportance();
idx = fm.setSelectedFeat(imp > 0);
[X, y] = SynEM.Util.getTrainingDataFrom([saveFolder, 'train'], labelType);
if exist('ensembleArgs','var') && ~isempty(ensembleArgs)
    classifier = SynEM.Classifier.BoostedEnsemble.train(X(:,idx), y, ...
        ensembleArgs{:});
else
    classifier = SynEM.Classifier.BoostedEnsemble.train(X(:,idx), y);
end
Util.save([saveFolder, 'classifier/ClassifierPred.mat'], classifier, imp);
end

