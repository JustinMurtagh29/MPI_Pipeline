function visualizeTestGP(filename, titleSuffix, filenameSuffix)

% This routine plots several metrics of a previous test result set.
%
% Use:
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se400.mat'], '400 Inducing Points, 26 Iterations', '.se400');
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se200.mat'], '200 Inducing Points, 20 Iterations', '.se200');
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se100-0.mat'], '100 Inducing Points, 20 Iterations', '.se100-0');
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se100-1.mat'], '100 Inducing Points, 40 Iterations', '.se100-1');

  load(filename);
  % global B; % This is expected to be a column vector of all predicted probabilites

  if ~isempty(titleSuffix)
    titleSuffix = [' (' titleSuffix ')'];
  end

  
  % %% PLOT true/false positive/negative curves per threshold
  % figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  
  % data = [];
  % for i=linspace(0,1, 101)
  %   [~, ~, ~, truePositives, trueNegatives, falsePositives, falseNegatives] = ...
  %     errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
  %   data = [data; [truePositives, trueNegatives, falsePositives, falseNegatives]];
  % end
  % plot(linspace(0,1, 101), data);
  % legend('truePositives', 'trueNegatives', 'falsePositives', 'falseNegatives');

  %% CALC precision/recall values per lower/upper threshold
  lower=[];
  for i=linspace(0,1, 101)
    [precision, recall, F1] = ...
      errorAnalysis((groundTruth.testLabels < 0), (predictedLabelProb <= i));
    lower = [lower; precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall), i]; % , sum(B(1:floor(i * 100 + 1)))
    % fprintf('Lower Threshold: %.2f   Precision: %f   Recall: %f   F1: %f   F0.5: %f\n', ...
    %   i, precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall));
  end
  upper=[];
  for i=linspace(0,1, 101)
    [precision, recall, F1] = ...
      errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
    upper = [upper; precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall), i]; % , sum(B(floor(i * 100 + 1):end))
    % fprintf('Upper Threshold: %.2f   Precision: %f   Recall: %f   F1: %f   F0.5: %f\n', ...
    %   i, precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall));
  end

  %% PLOT F1 score curves with marked maximum
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  plot(linspace(0,1,101), lower(:,3), 'b');
  plot(linspace(0,1,101), upper(:,3), 'r');
  legend({'lower','upper'},'Location','south');
  [m,idx] = max(lower(:,3));
  scatter((idx - 1)/100, m, 'x');
  text((idx - 1)/100 + 0.02, m + 0.02, ...
    sprintf('x: %.2f F: %.4f P: %.4f R: %.4f', (idx - 1)/100, m, lower(idx,1), lower(idx,2)));
  [m,idx] = max(upper(:,3));
  scatter((idx - 1)/100, m, 'x');
  text((idx - 1)/100 + 0.02, m + 0.02, ...
    sprintf('x: %.2f F: %.4f P: %.4f R: %.4f', (idx - 1)/100, m, upper(idx,1), upper(idx,2)));
  hold off
  title(['F1 curves' titleSuffix]);
  saveas(gcf, ['gp_f1' filenameSuffix '.pdf']);

  %% PLOT F0.5 score curves with marked maximum
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  plot(linspace(0,1,101), lower(:,4), 'b');
  plot(linspace(0,1,101), upper(:,4), 'r');
  legend({'lower','upper'},'Location','south');
  [m,idx] = max(lower(:,4));
  scatter((idx - 1)/100, m, 'x');
  text((idx - 1)/100 + 0.02, m + 0.02, ...
    sprintf('x: %.2f F: %.4f P: %.4f R: %.4f', (idx - 1)/100, m, lower(idx,1), lower(idx,2)));
  [m,idx] = max(upper(:,4));
  scatter((idx - 1)/100, m, 'x');
  text((idx - 1)/100 + 0.02, m + 0.02, ...
    sprintf('x: %.2f F: %.4f P: %.4f R: %.4f', (idx - 1)/100, m, upper(idx,1), upper(idx,2)));
  hold off
  title(['F0.5 curves' titleSuffix]);
  saveas(gcf, ['gp_f05' filenameSuffix '.pdf']);


  %% CALC thresholds with >= 0.95 precision and max F1
  lowerSubset = lower(find(lower(:,1) >= 0.95),:);
  [m,idx] = max(lowerSubset(:,3));
  lowerThreshold = lowerSubset(idx, 5);
  fprintf('%s Lower T: %.2f P: %.4f R: %.4f F1: %.4f\n', ...
    filenameSuffix, lowerSubset(idx, 5), lowerSubset(idx, 1), lowerSubset(idx, 2), lowerSubset(idx, 3));
  upperSubset = upper(find(upper(:,1) >= 0.97),:);
  [m,idx] = max(upperSubset(:,3));
  upperThreshold = upperSubset(idx, 5);
  fprintf('%s Upper T: %.2f P: %.4f R: %.4f F1: %.4f\n', ...
    filenameSuffix, upperSubset(idx, 5), upperSubset(idx, 1), upperSubset(idx, 2), upperSubset(idx, 3));

  
  %% PLOT precision/recall curve
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  plot(lower(1+1:end,2), lower(1+1:end,1), 'r');
  plot(upper(1:end-1,2), upper(1:end-1,1), 'b');
  legend({'lower', 'upper'},'Location','southwest');
  scatter(lower(floor(lowerThreshold * 100 + 1), 2), lower(floor(lowerThreshold * 100 + 1), 1), 'xr');
  scatter(upper(floor(upperThreshold * 100 + 1), 2), upper(floor(upperThreshold * 100 + 1), 1), 'xb');
  xlabel('recall');
  ylabel('precision');
  title(['Precision-Recall curves' titleSuffix]);
  hold off
  saveas(gcf, ['gp_precrec2' filenameSuffix '.pdf']);

  %% CALC nuisance/error
  data = [];
  for i=0:0.01:1
    for j=(i+0.01):0.01:1
      nuisance = sum(i < predictedLabelProb & predictedLabelProb < j) / length(predictedLabelProb);
      classSplit = [ones(sum(i <= predictedLabelProb), 1); zeros(sum(i > predictedLabelProb), 1)];
      splitError = sum(classSplit ~= (groundTruth.testLabels < 0)) / length(predictedLabelProb);
      classMerge = [zeros(sum(j < predictedLabelProb), 1); ones(sum(j >= predictedLabelProb), 1)];
      mergeError = sum(classMerge ~= (groundTruth.testLabels > 0)) / length(predictedLabelProb);
      data = [data; nuisance, splitError, mergeError];
    end
  end

  %% PLOT nuisance/error
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  scatter(data(:, 1), data(:, 2), 'rx');
  scatter(data(:, 1), data(:, 3), 'bx');
  hold off
  xlim([0 1]);
  ylim([0 1]);
  xlabel('Fraction of samples to be annotated');
  ylabel('Fraction of errors in test data');
  legend({'Split error', 'Merge error'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_error' filenameSuffix '.pdf']);

  %% CALC nuisance/error
  dataLower = [];
  for i=0:0.01:1
    samplesIdx = predictedLabelProb <= i;
    elimination = sum(samplesIdx) / length(predictedLabelProb);
    % samplesFullIdx = B <= i;
    % elimination = sum(samplesFullIdx) / length(B);
    FP_lower = sum(groundTruth.testLabels(samplesIdx) > 0) / sum(samplesIdx);
    dataLower = [dataLower; elimination, FP_lower];
  end
  dataUpper = [];
  for i=0:0.01:1
    samplesIdx = predictedLabelProb >= i;
    elimination = sum(samplesIdx) / length(predictedLabelProb);
    % samplesFullIdx = B >= i;
    % elimination = sum(samplesFullIdx) / length(B);
    FP_upper = sum(groundTruth.testLabels(samplesIdx) < 0) / sum(samplesIdx);
    dataUpper = [dataUpper; elimination, FP_upper];
  end

  %% PLOT nuisance/error
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  scatter(dataLower(:, 1), dataLower(:, 2), 'rx');
  scatter(dataUpper(:, 1), dataUpper(:, 2), 'bx');
  hold off
  xlim([0 1]);
  ylim([0 1]);
  xlabel('Fraction of test samples eliminated');
  ylabel('Fraction of errors in test data');
  legend({'Wrong splits/Missed merges', 'Wrong merges/Missed splits'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_error2' filenameSuffix '.pdf']);



  %% PLOT probabilities and confidences
  figure('PaperType', 'A4', 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'Renderer', 'painters', 'RendererMode', 'manual');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  paperSize = get(gcf,'PaperSize');
  set(gcf,'PaperPosition', [0, 0, paperSize(1), paperSize(2)]);
  [probTemp, idx] = sort(predictedLabelProb);
  hold on;
  idx1 = find(probTemp > lowerThreshold, 1, 'first');
  idx2 = find(probTemp > upperThreshold, 1, 'first');
  p = patch([0 idx1 idx1 0],[0 0 1 1], [1,0.67,0.62], 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], [0.46, 0.55, 1], 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  p = patch([idx2 length(probTemp) length(probTemp) idx2],[0 0 1 1], [0.5, 1, 0.42], 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  confidence(:,1) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence(:,2) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence = abs(confidence(:,1) - confidence(:,2));
  plot(confidence, 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  % plot(confidence(:,2), 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  tmp = median(binning(confidence, 100));
  plot(50+1:100:length(tmp)*100,tmp, 'Color', [0.6, 0.6, 0.6]);
  plot(probTemp, 'k');
  text(100, 0.02, num2str(idx1));
  text(idx1 + 100, 0.02, num2str(idx2 - idx1));
  text(idx2 + 100, 0.02, num2str(length(predictedLabelProb) - idx2));
  xlim([0 size(probTemp,1)]);
  xlabel('Edges sorted according to probability');
  ylabel('Probability according to classifier');
  legend({'Reject', 'Query', 'Accept', 'Uncertainty', 'Uncertainty Sliding Median (size=100)', 'Probability'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_stats' filenameSuffix '.pdf']);


  %% PLOT probabilities and confidences (query region only)
  figure('PaperType', 'A4', 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'Renderer', 'painters', 'RendererMode', 'manual');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  paperSize = get(gcf,'PaperSize');
  set(gcf,'PaperPosition', [0, 0, paperSize(1), paperSize(2)]);
  clear confidence;
  [probTemp, idx] = sort(predictedLabelProb);
  idx = find(lowerThreshold < probTemp & probTemp < upperThreshold);
  probTemp = probTemp(idx);
  hold on;
  confidence(:,1) = exp(feval(gpOptions.likfunc,[],ones(size(probTemp)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence(:,2) = exp(feval(gpOptions.likfunc,[],ones(size(probTemp)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence = abs(confidence(:,1) - confidence(:,2));
  plot(confidence, 'Color', 'b', 'LineWidth', 0.5);
  plot(probTemp, 'k', 'LineWidth', 2);
  posIdx = find(groundTruth.testLabels(idx) > 0);
  negIdx = find(groundTruth.testLabels(idx) < 0);
  scatter(posIdx, confidence(posIdx), 'bx');
  scatter(negIdx, confidence(negIdx), 'rx');
  xlim([0 size(probTemp,1)]);
  ylim([0 1]);
  xlabel('Edges sorted according to probability');
  ylabel('Probability according to classifier');
  legend({'Uncertainty', 'Probability', 'Label +1', 'Label -1'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_stats2' filenameSuffix '.pdf']);


  %% PLOT probabilities and actual frequencies
  figure('PaperType', 'A4', 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'Renderer', 'painters', 'RendererMode', 'manual');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  paperSize = get(gcf,'PaperSize');
  set(gcf,'PaperPosition', [0, 0, paperSize(1), paperSize(2)]);
  [probTemp, idx] = sort(predictedLabelProb);
  hold on;
  % confidence(:,1) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx)));
  % confidence(:,2) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx)));
  % plot(confidence(:,1), 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  % plot(confidence(:,2), 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  tmp = mean(binning(groundTruth.testLabels(idx) > 0, 100));
  bar((100/2)+1:100:length(tmp)*100,tmp);
  plot(probTemp, 'k');
  xlim([0 size(probTemp,1)]);
  xlabel('Edges sorted according to probability');
  ylabel('Probability according to classifier');
  legend({'Actual Frequency', 'Probability'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_stats3' filenameSuffix '.pdf']);

end


function fb = fBeta(beta, precision, recall)

  fb = (1 + beta*beta) * (precision * recall) / ((beta*beta * precision) + recall);

end

function bins = binning(A, binSize)

  if mod(length(A),binSize) ~= 0
    tmp = padarray(A, binSize - mod(length(A), binSize), 0, 'post');
  else
    tmp = A;
  end
  bins = reshape(tmp, [binSize ceil(length(A) / binSize)]);

end
