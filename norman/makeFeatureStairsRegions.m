% Plots for comparing the histograms of each feature in the 
% different training regions. One plot each for both labels.
%
% Use:
% load(pT.local(1,1,1).weightFile);
% load(pT.gp.initalGroundTruth);
% makeFeatureStairsRegions(regionData, weightNames)
function makeFeatureStairsRegions(regionData, weightNames)

  BIN_COUNT = 40;

  for i=[1:size(regionData(1).X, 2)]

    figure();
    set(gcf,'Visible', 'off');

    subplot(2,1,1);

    bins = {};
    centers = {};
    for j=1:length(regionData)
      positiveIndices = find(regionData(j).Y > 0);
      [bins{j}, centers{j}] = hist(regionData(j).X(positiveIndices,i), linspace(-1.5, 1.5, BIN_COUNT));
    end

    stairData = [];
    stairCenters = [];
    stairLabels = {};
    for j=1:length(bins)
      stairData = [stairData, (bins{j}/norm(bins{j}))'];
      stairCenters = [stairCenters, centers{j}'];
      stairLabels{j} = sprintf('Region %d', j);
    end

    stairs(stairCenters, stairData);
    legend(stairLabels);
    title(sprintf('Connected Feature %03d: %s', i, weightNames{i}));

    subplot(2,1,2);

    bins = {};
    centers = {};
    for j=1:length(regionData)
      negativeIndices = find(regionData(j).Y < 0);
      [bins{j}, centers{j}] = hist(regionData(j).X(negativeIndices,i), linspace(-1.5, 1.5, BIN_COUNT));
    end

    stairData = [];
    stairCenters = [];
    stairLabels = {};
    for j=1:length(bins)
      stairData = [stairData, (bins{j}/norm(bins{j}))'];
      stairCenters = [stairCenters, centers{j}'];
      stairLabels{j} = sprintf('Region %d', j);
    end

    stairs(stairCenters, stairData);
    legend(stairLabels);
    title(sprintf('Unconnected Feature %03d: %s', i, weightNames{i}));


    saveas(gcf, sprintf('feature%03d.pdf', i));
    disp(sprintf('feature%03d.pdf', i));
  end

end