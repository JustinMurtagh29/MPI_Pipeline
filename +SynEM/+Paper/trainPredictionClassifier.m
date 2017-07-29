function trainPredictionClassifier( featureFolder, ensembleArgs )
%TRAINPREDICTIONCLASSIFIER Train a classifier for prediction using only the
%necessary features.
% INPUT featureFolder: string
%           Path to folder where features are saved.
%           (see also Paper.SynEM.trainingPipeline saveFolder input)
%       ensembleArgs: (Optional) [2Nx1] cell
%           Cell array of name value pairs for
%           SynEM.Classifier.BoostedEnsemble.train.
% Author: Benedikt Staffler <benedikt.staffler@brain.mpg.de>

featureFolder = SynEM.Util.addFilesep(featureFolder);
m = load([featureFolder, 'classifier/LogitBoost.mat']);
imp = m.classifier.ens.predictorImportance();
[X, y, fm] = SynEM.Util.getTrainingDataFrom([featureFolder, 'train']);

%only use selected features
fm.setSelectedFeat(imp > 0);
idx = cell2mat(cellfun(@(x)x(:), fm.selectedFeat, 'uni', 0));
X = X(:, idx);

fprintf('[%s] Training classifier on %d features.\n', datestr(now), ...
    sum(idx));
if exist('ensembleArgs','var') && ~isempty(ensembleArgs)
    classifier = SynEM.Classifier.BoostedEnsemble.train(X, y, ...
        ensembleArgs{:});
else
    classifier = SynEM.Classifier.BoostedEnsemble.train(X, y);
end

fprintf('[%s] Saving classifier.\n', datestr(now));
save([featureFolder, 'classifier/LogitBoostPredFull.mat'], ...
    'classifier', '-v7.3');
classifier = compact(classifier);
classifier = classifier.calculatePredVar(); %#ok<NASGU>
save([featureFolder, 'classifier/LogitBoostPred.mat'], 'classifier');

end

