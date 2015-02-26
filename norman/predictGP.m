function [labelMean, labelCov, latentMean, latentCov, testLabelProb, post] = ...
  predictGP(testData, normFile, trainingData, trainingLabels, gpOptions)

  global GLOBAL_CODE_DIR;
  % Removing persistent values in @infFITC_EP
  clear infFITC_EP;

  if ~isempty(normFile)
    % normalize test data (variable weights)
    testData = normalizeDataForGP(testData, false, normFile);
  end

  initialTestLabels = ones(size(testData, 1), 1);

  if isreal(testData) == 0
    error('testData contains complex items.');
  end

  % gpml toolbox usage
  run([GLOBAL_CODE_DIR 'active/gpml/startup.m']);

  % gpOptions.hyp.cov(1:end-1) = gpOptions.hyp.cov(1:end-1) * 10;
  % Make predictions
  [labelMean, labelCov, latentMean, latentCov, lp, post] = ...
    gp(gpOptions.hyp, gpOptions.inffunc, gpOptions.meanfunc, gpOptions.covfunc, ...
      gpOptions.likfunc, trainingData, trainingLabels, testData, initialTestLabels);

  testLabelProb = exp(lp);
  
end

