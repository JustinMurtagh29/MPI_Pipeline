function runTestGP(parameter)

  groundTruth = load(parameter.gp.initalGroundTruth);
  gpOptions = load(parameter.gp.hyperParameter);

  outputFile = [parameter.gp.stateFolder 'testResults.' datestr(clock, 30) '.mat'];
  disp(['Output to: ' outputFile]);

  [labelMean, labelCov, latentMean, latentCov, predictedLabelProb, post] = ...
    predictGP(groundTruth.testData, parameter.gp.normValues, groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);
  % Small version
  % predictGP(groundTruth.testData, parameter.gp.normValues, groundTruth.trainingData(1:1000,:), groundTruth.trainingLabels(1:1000,:), gpOptions);
  
  save(outputFile, 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb', 'post');

  prec_rec(predictedLabelProb, groundTruth.testLabels)

  % data = [];
  % lower=[];
  % for i=linspace(0,1,100)
  %   [precision, recall, F1, truePositives, trueNegatives, falsePositives, falseNegatives] = ...
  %     errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
  %   lower = [lower; precision, recall];
  %   data = [data; [truePositives, trueNegatives, falsePositives, falseNegatives]];
  %   disp(sprintf('Threshold: %.2f   Precision: %f   Recall: %f   F1: %f', i, precision, recall, F1));
  % end
  % plot(data);
  % legend('truePositives', 'trueNegatives', 'falsePositives', 'falseNegatives');


  % upper=[];
  % for i=linspace(0,1,100)
  %   [precision, recall, F1, truePositives, trueNegatives, falsePositives, falseNegatives] = ...
  %     errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
  %   upper = [upper; precision, recall];
  %   disp(sprintf('Threshold: %.2f   Precision: %f   Recall: %f   F1: %f', i, precision, recall, F1));
  % end

  % hold on
  % scatter(lower(:,2), lower(:,1), [],linspace(0,1,size(lower, 1)));
  % scatter(upper(:,2), upper(:,1), [],linspace(1,0,size(upper, 1)));
  % hold off

  % plot(linspace(0,1,100), [lower upper]);


  % % Plot probabilities and covariance
  % plot_variance = @(x,lower,upper,color) fill([x;flip(x)],[upper;flip(lower)],color,'EdgeColor',color);
  % plot_variance([1:7012]', data(:,1) + abs(data(:,2) / 2), data(:,1) - abs(data(:,2) / 2), [0.8 0.8 0.8]);
  % hold on
  % plot([1:7012],data(:,1));
  % hold off

end