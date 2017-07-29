function trainingPipeline( dataFolder, saveFolder, fm, pTest, options, ...
    cluster )
%TRAININGPIPELINE Training pipeline for SynEM classifiers.
% INPUT dataFolder: string
%           Path to training data folder.
%       saveFolder: string
%           Path to folder where features are saved.
%       featureMap: Interface.FeatureMap object
%       pTest: (Optional) struct
%           Segmentation parameter struct used to load the information for
%           cube 67 (test cube).
%       options: (Optional) Struct
%           Struct with further options
%           'labelType': string
%               See 'type' in SynEM.Util.getTrainingDataFrom.
%               (Default: 'direction')
%           'discardUndecided': logical
%               Discard interfaces marked as undecided during training.
%               see also SynEM.Training.calculateFeaturesForTrainingData
%               (Default: false)
%           'ensembleArgs: [2Nx1] cell
%               Cell array of name value pairs for
%               SynEM.Classifier.BoostedEnsemble.train.
%           'skipFeatCalc': logical
%               Skip the feature calculation and only use the
%               features found in the 'test' and 'val' subfolders.
%               (Default: false)
%           'className': string
%               Filename for the classifier output file (without .mat
%               ending).
%               (Default: LogitBoost)
%           'trainPredClassifier': logical
%               Flag indicating that also a classifier only on the relevant
%               features should be trained.
%               (Default: true)
%           'overwriteInputs': logical
%               Overwrite the TrainingPipelineInputs.mat in the training
%               main folder if it already exists.
%               (Default: false)
%           'loadInputParameters': logical
%               Flag to load the TrainingPipelineInputs.mat from a previous
%               run. If true the only the saveFolder has to be specified.
%               Everything except the options are then replaced by the
%               content of 'TrainingPipelineInputs.mat' from the
%               saveFolder.
%               (Default: false)
%       cluster: (Optional) parallel.cluster object
%           The cluster object used for parallel computation.
%           (Default: getCluster('cpu') or local computation if the former
%           does not exist)
% see also SynEM.Training.calculateFeaturesForTrainingData
%
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

%options
if ~exist('options', 'var') || isempty(options)
    options = struct;
end
options = defaultOptions(options);

%load previous input parameters if specified
if options.loadInputParameters
    m = load(fullfile(fullfile(saveFolder, 'TrainingPipelineInputs.mat')));
    dataFolder = m.dataFolder;
    fm = m.fm;
    fm.setSelectedFeat(); %set to all features
    pTest = m.pTest;
end

dataFolder = SynEM.Util.addFilesep(dataFolder);
saveFolder = SynEM.Util.addFilesep(saveFolder);

%store git hash
try
    git_hash = Util.getGitHash();
catch
    git_hash = [];
    warning('Unable to get git commit hash.');
end

%store all inputs for reproducability
if ~exist(fullfile(saveFolder, 'TrainingPipelineInputs.mat'), 'file') || ...
        options.overwriteInputs
    Util.log('Pipeline run parameters are stored in %s.', ...
        fullfile(saveFolder, 'TrainingPipelineInputs.mat'));
    Util.save(fullfile(saveFolder, 'TrainingPipelineInputs.mat'), ...
        dataFolder, saveFolder, fm, pTest, git_hash, options);
else
    Util.log(['Training pipeline parameters from previous run were' ...
        ' found and will not be overwritten.']);
end

if ~exist(saveFolder,'dir')
    mkdir(saveFolder);
end
Util.log('Results are stored in %s.', saveFolder);

% feature calculation
if ~options.skipFeatCalc
    Util.log('Calculating features.');
    if exist('cluster', 'var') && ~isempty(cluster)
        job = SynEM.Training.calculateFeaturesForTrainingData( ...
            dataFolder, saveFolder, fm, true, cluster, ...
            options.discardUndecided);
    else
        job = SynEM.Training.calculateFeaturesForTrainingData( ...
            dataFolder, saveFolder, fm, true, [], options.undecidedList);
    end

    if ~isempty(job) % if calculation was done on cluster
        Util.log('Waiting for job %d results.', job.Id);
        wait(job);
        errIdx = Cluster.getIdxOfTasksWithError(job);
        if ~isempty(errIdx)
            error(['The feature calculation job with id %d has errors.' ...
                ' Aborting training pipeline ...'], job.Id);
        end
    end
end
SynEM.Util.groupTrainingData(saveFolder);

% classifier training
Util.log('Training classifier.');
[X, y] = SynEM.Util.getTrainingDataFrom([saveFolder, 'train'], ...
    options.labelType);

if ~isempty(options.ensembleArgs)
    classifier = SynEM.Classifier.BoostedEnsemble.train(X, y, ...
        options.ensembleArgs{:});
else
    classifier = SynEM.Classifier.BoostedEnsemble.train(X, y);
end
classifier.options.fm = fm;
classifier.options.git_hash = git_hash;

%classifier saving
Util.log('Saving classifier to %s.', [saveFolder, 'classifier']);
if ~exist([saveFolder, 'classifier'], 'dir')
    mkdir([saveFolder, 'classifier']);
end
Util.save(fullfile(saveFolder, 'classifier', ...
    [options.className 'Full.mat']), ...
    classifier);
classifier = compact(classifier);
classifier = classifier.calculatePredVar();
Util.save(fullfile(saveFolder, 'classifier', ...
    [options.className '.mat']), ...
    classifier);

%classifier evaluation
Util.log('Evaluating validation set performance.');

%training performance
[~, scoresTrain] = classifier.predict(X);
[rpTrain, ~, aucTrain] = SynEM.Eval.interfaceRP(y, scoresTrain);

%validation performance
if strcmp(options.labelType, 'augment')
    %use single lables if trained using augmented labels
    options.labelType = 'single';
end
[X, y] = SynEM.Util.getTrainingDataFrom([saveFolder, 'val'], ...
    options.labelType);
[~, scoresVal] = classifier.predict(X);
[rpVal, ~, aucVal] = SynEM.Eval.interfaceRP(y, scoresVal);

if ~exist('pTest', 'var') || isempty(pTest)
    Util.save(fullfile(saveFolder, 'classifier', ...
        ['result_' options.className '.mat']), ...
        scoresVal, rpVal);
else
    Util.log('Evaluating test set performance.');
    scoresTest = SynEM.Seg.predictCube(pTest, 67, fm, classifier);
    scoresTest = reshape(scoresTest, [], 2);

    savePath = fullfile(saveFolder, 'classifier', ...
        ['result_' options.className '.mat']);
    try %to load the test set on gaba
        testSet = SynEM.Test.loadInterfaceGT( ...
            '/u/bstaffle/data/SynEM/data', true, true);
        group = testSet.group;
        toKeep = testSet.toKeep;
        [rpTest, ~, aucTest] = SynEM.Eval.interfaceRP(group, ...
            max(scoresTest(toKeep,:), [], 2));
        Util.save(savePath, scoresTest, rpTest, scoresVal, rpVal, ...
            scoresTrain, rpTrain, aucTrain, aucVal, aucTest);
    catch err
        warning('Could not load the test set data due to %s.', ...
            err.message);
        Util.save(savePath, scoresTest, scoresVal, rpVal, scoresTrain, ...
            rpTrain, aucTrain, aucVal);
    end
    Util.log('Classifier results saved to %s.', savePath);
end

%training prediction classifier
if options.trainPredClassifier
    Util.log('Training prediction classifier.');
    imp = classifier.ens.predictorImportance();
    idx = fm.setSelectedFeat(imp > 0);
    if ~isempty(options.ensembleArgs)
        classifierPred = SynEM.Classifier.BoostedEnsemble.train( ...
            X(:,idx), y, options.ensembleArgs{:});
    else
        classifierPred = SynEM.Classifier.BoostedEnsemble.train( ...
            X(:,idx), y);
    end
    classifierPred.options.fm = fm;
    classifierPred.options.imp = imp;
    classifierPred.options.git_hash = git_hash;

    classifierPred = classifierPred.compact();

    savePath = fullfile(saveFolder, 'classifier', ...
        [options.className 'Pred.mat']);
    Util.log('Saving prediction classifiert to %s.', savePath);
    Util.save(savePath, classifierPred);
end

end

function options = defaultOptions(options)

if ~isfield(options, 'labelType')
    options.labelType = 'direction';
end
if ~isfield(options, 'discardUndecided')
    options.discardUndecided = false;
end
if ~isfield(options, 'ensembleArgs')
    options.ensembleArgs = [];
end
if ~isfield(options, 'overwriteInputs')
    options.overwriteInputs = false;
end
if ~isfield(options, 'loadInputParameters')
    options.loadInputParameters = false;
end
if ~isfield(options, 'skipFeatCalc')
    options.skipFeatCalc = false;
end
if ~isfield(options, 'className')
    options.className = 'LogitBoost';
else
    [~,options.className] = fileparts(options.className);
end
if ~isfield(options, 'trainPredClassifier')
    options.trainPredClassifier = true;
end
end
