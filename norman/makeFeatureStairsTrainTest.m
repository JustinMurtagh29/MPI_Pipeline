function makeFeatureStairsTrainTest(trainingData, trainingLabels, testData, testLabels, weightNames, hyp)

% Plots for comparing the histograms of each feature in training and test set.
%
% Use:
% load(pT.local(1,1,1).weightFile);
% load(pT.gp.initalGroundTruth);
% load(pT.gp.hyperParameter);
% makeFeatureStairsTrainTest(trainingData, trainingLabels, testData, testLabels, weightNames, hyp)

  BIN_COUNT = 40;

  hyp = hyp.cov(1:end-1);
  [hyp, featureIdx] = sort(hyp);

  for i=[1:size(trainingData, 2)]

    j = featureIdx(i);

    stairLabels = {};
    figure();
    set(gcf,'Visible', 'off');

    % training +1
    idx = find(trainingLabels > 0);
    [bins, centers] = hist(trainingData(idx,i), linspace(-1.5, 1.5, BIN_COUNT));
    
    stairs(centers', (bins/norm(bins))', '-r');
    hold on;
    stairLabels{1} = sprintf('Training +1');


    % training -1
    idx = find(trainingLabels < 0);
    [bins, centers] = hist(trainingData(idx,i), linspace(-1.5, 1.5, BIN_COUNT));
    
    stairs(centers', (bins/norm(bins))', '-b');
    hold on;
    stairLabels{2} = sprintf('Training -1');

    % test +1
    idx = find(testLabels > 0);
    [bins, centers] = hist(testData(idx,i), linspace(-1.5, 1.5, BIN_COUNT));
    
    stairs(centers', (bins/norm(bins))', ':r');
    hold on;
    stairLabels{3} = sprintf('Test +1');


    % test -1
    idx = find(testLabels < 0);
    [bins, centers] = hist(testData(idx,i), linspace(-1.5, 1.5, BIN_COUNT));
    
    stairs(centers', (bins/norm(bins))', ':b');
    hold on;
    stairLabels{4} = sprintf('Test -1');


    legend(stairLabels);
    title(sprintf('Feature %03d: %s hyp=%4.5f', i, weightNames{i}, hyp(j)));

    saveas(gcf, sprintf('feature%03d.pdf', j));
    disp(sprintf('feature%03d.pdf %i/%i', j, i, size(trainingData, 2)));
  end

end