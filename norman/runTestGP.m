function runTestGP(parameter)

  groundTruth = load(parameter.gp.initalGroundTruth);
  gpOptions = load(parameter.gp.hyperParameter);

  outputFile = [parameter.gp.stateFolder 'testResults.' datestr(clock, 30) '.mat'];
  disp(['Output to: ' outputFile]);

  [labelMean labelCov latentMean latentCov predictedLabelProb post] = ...
    predictGP(groundTruth.testData, parameter.gp.normValues, groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);
  % Small version
  % predictGP(groundTruth.testData, parameter.gp.normValues, groundTruth.trainingData(1:1000,:), groundTruth.trainingLabels(1:1000,:), gpOptions);
  
  save(outputFile, 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb', 'post');

  for i=[0:.02:1]
    [precision, recall, F1] = errorAnalysis(groundTruth.testLabels, predictedLabelProb >= i);
    disp(sprintf('Threshold: %f   Precision: %f   Recall: %f   F1: %f', i, precision, recall, F1));
  end
end