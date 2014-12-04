function runTestGP(parameter, hyp)

% This routine predicts classification for the test set.
% It stores the results for further analysis.
% 
% Use:
% a) Default
%    runTestGP(pT);
% b) Use a different hyperParameter set
%    runTestGP(pT, hyp);

  groundTruth = load(parameter.gp.initalGroundTruth);
  gpOptions = load(parameter.gp.hyperParameter);
  
  if nargin == 2
    gpOptions.hyp = hyp;
  end

  outputFile = [parameter.gp.stateFolder 'testSet.' datestr(clock, 30) '.mat'];
  disp(['Output to: ' outputFile]);

  %% RUN prediction
  [labelMean, labelCov, latentMean, latentCov, predictedLabelProb] = ...
    predictGP(groundTruth.testData, [], groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);
  
  save(outputFile, 'groundTruth', 'gpOptions', 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb');

end
