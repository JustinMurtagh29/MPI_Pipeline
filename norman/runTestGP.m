function runTestGP(parameter)

  groundTruth = load(parameter.gp.initalGroundTruth);
  gpOptions = load(parameter.gp.hyperParameter);

  [labelMean labelCov latentMean latentCov predictedLabelProb] = ...
    testGP(groundTruth.testData, parameter.gp.normValues, groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);
  
  save([parameter.gp.saveFolder 'testResults.mat'], 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb');

  errorAnalysis(groundTruth.testLabels, predictedLabelProb > 0.5)

end