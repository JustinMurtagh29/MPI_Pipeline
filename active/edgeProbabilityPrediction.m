function edgeProbabilityPrediction(weightFile, normFile, groundTruthFile, hyperFile, probFile)

% load data
load(weightFile); % test data to predict
load(groundTruthFile); % training data generated with prepareTrainingData
load(hyperFile); % all parameter for GP (e.g. hyp, meanfunc)

% normalize test data (variable weights)
weights = normalizeDataForGP(weights, false, normFile);

% gpml toolbox usage
run('./active/gpml/startup.m');
% Make predictions
[labelMean labelCov latentMean latentCov lp post] = gp(hyp, inffunc, meanfunc, covfunc, likfunc, trainingData, trainingLabels, weights, ones(size(weights,1), 1));

% save preditions to designated file
prob = exp(lp);
save(probFile, 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'prob', 'post');

end

