function runTestActiveGP(parameter, testSetFilename)

% Same as runTestGP, but repartitions training and test
% prior to prediction.

  load(testSetFilename);

  outputFile = [parameter.gp.stateFolder 'testSetActive.' datestr(clock, 30) '.mat'];
  disp(['Output to: ' outputFile]);

  %% Repartition training/test data
  lowerThreshold = 0.21;
  upperThreshold = 0.90;

  % Find test samples in query region
  idx = find(lowerThreshold <= predictedLabelProb & predictedLabelProb <= upperThreshold);

  % % Calc confidences (size of confidence interval 95%)
  % confidences = abs(...
  %   exp(feval(gpOptions.likfunc,[],ones(size(idx)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx))) - ...
  %   exp(feval(gpOptions.likfunc,[],ones(size(idx)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx))));

  % % Sort by confidences
  % [~, idx2] = sort(confidences);
  
  % % Pick the `k` most unconfident samples
  % k = 100;
  % idx = idx(idx2(end-k:end));

  idx2 = randperm(length(idx), 100);
  idx = idx(idx2);

  % Repartition
  groundTruth.trainingData = [groundTruth.trainingData; groundTruth.testData(idx,:)];
  groundTruth.trainingLabels = [groundTruth.trainingLabels; groundTruth.testLabels(idx,:)];

  groundTruth.testData = groundTruth.testData(setdiff(1:end,idx),:);
  groundTruth.testLabels = groundTruth.testLabels(setdiff(1:end,idx),:);


  %% RUN prediction
  [labelMean, labelCov, latentMean, latentCov, predictedLabelProb] = ...
    predictGP(groundTruth.testData, [], groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);

  save(outputFile, 'groundTruth', 'gpOptions', 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb');

end