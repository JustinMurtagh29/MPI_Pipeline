function [precision, recall, F1] = errorAnalysis(testLabels, predictedLabels)

  if size(testLabels) ~= size(predictedLabels)
    error('testLabels and predictedLabels don''t line up.');
  end

  % Set testLabels to either 0 or 1 (not -1 and 1)
  testLabels = testLabels > 0;

  positives = predictedLabels == 1;
  negatives = predictedLabels == 0;

  trues = positives == testLabels;
  falses = positives ~= testLabels;

  truePositives = sum(trues & positives);
  trueNegatives = sum(trues & negatives);
  falsePositives = sum(falses & positives);
  falseNegatives = sum(falses & negatives);

  if truePositives + trueNegatives + falsePositives + falseNegatives ~= size(predictedLabels, 1)
    error('Some samples have not been covered.')
  end

  % Precision
  if truePositives + falsePositives > 0
    precision = truePositives / (truePositives + falsePositives);
  else
    precision = 0;
  end

  % Recall
  if truePositives + falseNegatives > 0
    recall = truePositives / (truePositives + falseNegatives);
  else
    recall = 0;
  end

  % F1 Score
  if precision + recall > 0
    F1 = (2 * precision * recall) / (precision + recall);
  else
    F1 = 0;
  end

end