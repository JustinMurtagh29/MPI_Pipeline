function [labelMean labelCov latentMean latentCov testLabelProb] = ...
  predictGP(testData, normFile, trainingData, trainingLabels, gpOptions)

  % normalize test data (variable weights)
  testData = normalizeDataForGP(testData, false, normFile);
  initialTestLabels = ones(size(testData, 1), 1);

  % gpml toolbox usage
  run('/zdata/manuel/code/active/gpml/startup.m');

  % Make predictions
  [labelMean labelCov latentMean latentCov lp] = ...
    gp(gpOptions.hyp, gpOptions.inffunc, gpOptions.meanfunc, gpOptions.covfunc, ...
      gpOptions.likfunc, trainingData, trainingLabels, testData, initialTestLabels);

  testLabelProb = exp(lp);
  
end

