function [precision, recall, F1] = errorAnalysis(testLabels, predictedLabels)

  positives = predictedLabels == 1;
  negatives = predictedLabels == 0;

  trues = positives == predictedLabels;
  falses = positives ~= predictedLabels;

  truePositives = sum(trues & positives);
  trueNegatives = sum(trues & negatives);
  falsePositives = sum(falses & positives);
  falseNegatives = sum(falses & negatives);

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