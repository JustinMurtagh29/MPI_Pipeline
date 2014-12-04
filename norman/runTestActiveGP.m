function runTestActiveGP(parameter)

% Same as runTestGP, but repartitions training and test
% prior to prediction.

  groundTruth = load(parameter.gp.initalGroundTruth);
  gpOptions = load(parameter.gp.hyperParameter);

  outputFile = [parameter.gp.stateFolder 'testSetActive.' datestr(clock, 30) '.mat'];
  disp(['Output to: ' outputFile]);


  %% Repartition training/test data
  idx = randperm(floor(0.5 * size(groundTruth.testData, 1)));

  groundTruth.trainingData = [groundTruth.trainingData; groundTruth.testData(idx,:)];
  groundTruth.trainingLabels = [groundTruth.trainingLabels; groundTruth.testLabels(idx,:)];

  groundTruth.testData = groundTruth.testData(setdiff(1:end,idx),:);
  groundTruth.testLabels = groundTruth.testLabels(setdiff(1:end,idx),:);

  %% RUN prediction
  [labelMean, labelCov, latentMean, latentCov, predictedLabelProb] = ...
    predictGP(groundTruth.testData, [], groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);

  save(outputFile, 'groundTruth', 'gpOptions', 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb');

end