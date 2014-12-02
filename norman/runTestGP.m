function runTestGP(parameter, hyp)

  groundTruth = load(parameter.gp.initalGroundTruth);
  gpOptions = load(parameter.gp.hyperParameter);
  
  if nargin == 2
      gpOptions.hyp = hyp;
  end

  outputFile = [parameter.gp.stateFolder 'testResults.' datestr(clock, 30) '.mat'];
  disp(['Output to: ' outputFile]);

  [labelMean, labelCov, latentMean, latentCov, predictedLabelProb] = ...
    predictGP(groundTruth.testData, parameter.gp.normValues, groundTruth.trainingData, groundTruth.trainingLabels, gpOptions);
  
  % save(outputFile, 'labelMean', 'labelCov', 'latentMean', 'latentCov', 'predictedLabelProb');

  %% PLOT precision/recall curve
  prec_rec(predictedLabelProb, groundTruth.testLabels)

  
  %% PLOT true/false positive/negative curves per threshold
  figure;
  
  data = [];
  for i=linspace(0,1, 101)
    [~, ~, ~, truePositives, trueNegatives, falsePositives, falseNegatives] = ...
      errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
    data = [data; [truePositives, trueNegatives, falsePositives, falseNegatives]];
  end
  plot(linspace(0,1, 101), data);
  legend('truePositives', 'trueNegatives', 'falsePositives', 'falseNegatives');

  %% CALC precision/recall values per lower/upper threshold
  lower=[];
  for i=linspace(0,1, 101)
    [precision, recall, F1] = ...
      errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
    lower = [lower; precision, recall, F1];
    disp(sprintf('Threshold: %.2f   Precision: %f   Recall: %f   F1: %f', i, precision, recall, F1));
  end
  upper=[];
  for i=linspace(0,1, 101)
    [precision, recall, F1] = ...
      errorAnalysis((groundTruth.testLabels < 0), (predictedLabelProb <= i));
    upper = [upper; precision, recall, F1];
    disp(sprintf('Threshold: %.2f   Precision: %f   Recall: %f   F1: %f', i, precision, recall, F1));
  end

  %% PLOT F1 score curves with marked maximum
  figure;
  hold on
  plot(linspace(0,1,101), upper(:,3), 'r');
  plot(linspace(0,1,101), lower(:,3), 'b');
  legend({'upper','lower'},'Location','south');
  [m,idx] = max(upper(:,3));
  scatter((idx - 1)/100, m, 'x');
  [m,idx] = max(lower(:,3));
  scatter((idx - 1)/100, m, 'x');
  hold off

  
  %% PLOT precision/recall curve
  figure;
  hold on
  plot(lower(:,2), lower(:,1), 'r');
  plot(upper(:,2), upper(:,1), 'b');
  hold off
  legend({'lower', 'upper'},'Location','southwest');
  xlabel('recall');
  ylabel('precision');
  title('Precision-Recall curves');


  %% PLOT precision/recall scattered curves
  figure;
  hold on
  scatter(lower(:,2), lower(:,1), [], linspace(0,1,size(lower, 1)), 'x');
  scatter(upper(:,2), upper(:,1), [], linspace(1,0,size(upper, 1)), 'o');
  legend({'lower', 'upper'},'Location','southwest');
  hold off
  colorbar;
  title('Precision-Recall scattered curves');


  %% PLOT probabilities and confidences
  [probTemp, idx] = sort(predictedLabelProb);
  figure('Units', 'centimeters', 'Position', [0 0 29.7 21], ...
  'PaperType', 'A4', 'PaperUnits', 'centimeters', 'PaperPosition', [0 0 29.7 21], 'PaperOrientation', 'portrait');
  plot(probTemp, 'k', 'LineWidth', 2);
  hold on;
  idx1 = find(probTemp > .38, 1, 'first');
  idx2 = find(probTemp > .52, 1, 'first');
  p = patch([0 idx1 idx1 0],[0 0 1 1], 'r', 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], 'b', 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  p = patch([idx2 length(probTemp) length(probTemp) idx2],[0 0 1 1], 'g', 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  x = 1:length(probTemp);
  confidence(:,1) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence(:,2) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx)));
  plot(confidence(:,1), 'c');
  plot(confidence(:,2), 'c');
  xlim([0 size(probTemp,1)]);
  xlabel('Edges sorted according to probability');
  ylabel('Probability according to classifier');
  legend({'Probability', 'Reject', 'Query', 'Accept', 'Confidence'}, 'Location', 'Northwest');
  title('Classification statistics');

end