function visualizeTestGP(filename, titleSuffix, filenameSuffix)

% This routine plots several metrics of a previous test result set.
%
% Use:
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se400.mat'], '400 Inducing Points, 26 Iterations', '.se400');
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se200.mat'], '200 Inducing Points, 20 Iterations', '.se200');
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se100-0.mat'], '100 Inducing Points, 20 Iterations', '.se100-0');
% close all; visualizeTestGP([pT.gp.stateFolder 'testSet.se100-1.mat'], '100 Inducing Points, 40 Iterations', '.se100-1');

  load(filename);

  if ~isempty(titleSuffix)
    titleSuffix = [' (' titleSuffix ')'];
  end


  % %% PLOT precision/recall curve
  % prec_rec(predictedLabelProb, groundTruth.testLabels);

  
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
    lower = [lower; precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall)];
    fprintf('Lower Threshold: %.2f   Precision: %f   Recall: %f   F1: %f   F0.5: %f\n', ...
      i, precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall));
  end
  upper=[];
  for i=linspace(0,1, 101)
    [precision, recall, F1] = ...
      errorAnalysis((groundTruth.testLabels > 0), (predictedLabelProb >= i));
    upper = [upper; precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall)];
    fprintf('Upper Threshold: %.2f   Precision: %f   Recall: %f   F1: %f   F0.5: %f\n', ...
      i, precision, recall, fBeta(1, precision, recall), fBeta(0.5, precision, recall));
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
  lowerCutoff = (idx - 1)/100;
  scatter((idx - 1)/100, m, 'x');
  text((idx - 1)/100 + 0.02, m + 0.02, ...
    sprintf('x: %.2f F: %.4f P: %.4f R: %.4f', (idx - 1)/100, m, lower(idx,1), lower(idx,2)));
  [m,idx] = max(upper(:,4));
  upperCutoff = (idx - 1)/100;
  scatter((idx - 1)/100, m, 'x');
  text((idx - 1)/100 + 0.02, m + 0.02, ...
    sprintf('x: %.2f F: %.4f P: %.4f R: %.4f', (idx - 1)/100, m, upper(idx,1), upper(idx,2)));
  hold off
  title(['F0.5 curves' titleSuffix]);
  saveas(gcf, ['gp_f05' filenameSuffix '.pdf']);

  
  %% PLOT precision/recall curve
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  plot(lower(:,2), lower(:,1), 'r');
  plot(upper(:,2), upper(:,1), 'b');
  hold off
  legend({'lower', 'upper'},'Location','southwest');
  xlabel('recall');
  ylabel('precision');
  title(['Precision-Recall curves' titleSuffix]);
  saveas(gcf, ['gp_precrec2' filenameSuffix '.pdf']);


  %% PLOT precision/recall curves
  figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  hold on
  plot(linspace(0,1,101), lower(:,1), 'r');
  plot(linspace(0,1,101), lower(:,2), 'b');
  plot(linspace(0,1,101), upper(:,1), 'r--');
  plot(linspace(0,1,101), upper(:,2), 'b--');
  hold off
  legend({'Precision -1', 'Recall -1', 'Precision +1', 'Recall +1'},'Location','southwest');
  xlabel('recall');
  ylabel('precision');
  title(['Precision-Recall curves' titleSuffix]);
  saveas(gcf, ['gp_precrec' filenameSuffix '.pdf']);


  % %% PLOT precision/recall scattered curves
  % figure('PaperType', 'A4', 'PaperOrientation', 'portrait', 'PaperPositionMode', 'auto');
  % hold on
  % scatter(lower(:,2), lower(:,1), [], linspace(0,1,size(lower, 1)), 'x');
  % scatter(upper(:,2), upper(:,1), [], linspace(1,0,size(upper, 1)), 'o');
  % legend({'lower', 'upper'},'Location','southwest');
  % hold off
  % colorbar;
  % title(['Precision-Recall scattered curves' titleSuffix]);


  %% PLOT probabilities and confidences
  figure('PaperType', 'A4', 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'Renderer', 'painters', 'RendererMode', 'manual');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  paperSize = get(gcf,'PaperSize');
  set(gcf,'PaperPosition', [0, 0, paperSize(1), paperSize(2)]);
  [probTemp, idx] = sort(predictedLabelProb);
  hold on;
  idx1 = find(probTemp > lowerCutoff, 1, 'first');
  idx2 = find(probTemp > upperCutoff, 1, 'first');
  p = patch([0 idx1 idx1 0],[0 0 1 1], [1,0.67,0.62], 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], [0.46, 0.55, 1], 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  p = patch([idx2 length(probTemp) length(probTemp) idx2],[0 0 1 1], [0.5, 1, 0.42], 'EdgeColor', 'none');
  set(p,'FaceAlpha',0.5);
  x = 1:length(probTemp);
  confidence(:,1) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence(:,2) = exp(feval(gpOptions.likfunc,[],ones(size(latentMean)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx)));
  plot(abs(confidence(:,1) - confidence(:,2)), 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  % plot(confidence(:,2), 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  plot(probTemp, 'k');
  xlim([0 size(probTemp,1)]);
  xlabel('Edges sorted according to probability');
  ylabel('Probability according to classifier');
  legend({'Reject', 'Query', 'Accept', 'Uncertainty', 'Probability'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_stats' filenameSuffix '.pdf']);

  %% PLOT probabilities and confidences
  figure('PaperType', 'A4', 'PaperPositionMode', 'auto', 'PaperUnits', 'centimeters', ...
    'Renderer', 'painters', 'RendererMode', 'manual');
  set(gcf,'PaperSize',fliplr(get(gcf,'PaperSize')));
  paperSize = get(gcf,'PaperSize');
  set(gcf,'PaperPosition', [0, 0, paperSize(1), paperSize(2)]);
  clear confidence;
  [probTemp, idx] = sort(predictedLabelProb);
  idx = find(lowerCutoff < probTemp & probTemp < upperCutoff);
  probTemp = probTemp(idx);
  hold on;
  % idx1 = find(probTemp > lowerCutoff, 1, 'first');
  % idx2 = find(probTemp > upperCutoff, 1, 'first');
  % p = patch([0 idx1 idx1 0],[0 0 1 1], [1,0.67,0.62], 'EdgeColor', 'none');
  % set(p,'FaceAlpha',0.5);
  % p = patch([idx1 idx2 idx2 idx1],[0 0 1 1], [0.46, 0.55, 1], 'EdgeColor', 'none');
  % set(p,'FaceAlpha',0.5);
  % p = patch([idx2 length(probTemp) length(probTemp) idx2],[0 0 1 1], [0.5, 1, 0.42], 'EdgeColor', 'none');
  % set(p,'FaceAlpha',0.5);
  confidence(:,1) = exp(feval(gpOptions.likfunc,[],ones(size(probTemp)),latentMean(idx)-2*sqrt(latentCov(idx)), latentCov(idx)));
  confidence(:,2) = exp(feval(gpOptions.likfunc,[],ones(size(probTemp)),latentMean(idx)+2*sqrt(latentCov(idx)), latentCov(idx)));
  plot(abs(confidence(:,1) - confidence(:,2)), 'Color', 'b', 'LineWidth', 0.5);
  % plot(confidence(:,2), 'Color', [0.95, 0.95, 0.95], 'LineWidth', 0.5);
  plot(probTemp, 'k');
  xlim([0 size(probTemp,1)]);
  ylim([0 1]);
  xlabel('Edges sorted according to probability');
  ylabel('Probability according to classifier');
  legend({'Uncertainty', 'Probability'}, 'Location', 'Northwest');
  title(['Classification statistics' titleSuffix]);
  saveas(gcf, ['gp_stats2' filenameSuffix '.pdf']);

end


function fb = fBeta(beta, precision, recall)

  fb = (1 + beta*beta) * (precision * recall) / ((beta*beta * precision) + recall);

end
