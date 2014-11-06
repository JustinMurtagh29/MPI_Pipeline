function makeFeatureStairs(regionData, weightNames)

  BIN_COUNT = 40;

  parfor i=[1:size(regionData(1).X, 2)]

    figure();
    set(gcf,'Visible', 'off');

    subplot(2,1,1);

    bins = {};
    for j=1:length(regionData)
      positiveIndices = find(regionData(j).Y > 0);
      bins{j} = hist(regionData(j).X(positiveIndices,i), BIN_COUNT);
    end

    stairData = [];
    stairLabels = {};
    for j=1:length(bins)
      stairData = [stairData, (bins{j}/norm(bins{j}))'];
      stairLabels{j} = sprintf('Region %d', j);
    end

    stairs([1:BIN_COUNT], stairData);
    legend(stairLabels);
    title(sprintf('Connected Feature %03d: %s', i, weightNames{i}));

    subplot(2,1,2);

    bins = {};
    for j=1:length(regionData)
      negativeIndices = find(regionData(j).Y < 0);
      bins{j} = hist(regionData(j).X(negativeIndices,i), BIN_COUNT);
    end

    stairData = [];
    stairLabels = {};
    for j=1:length(bins)
      stairData = [stairData, (bins{j}/norm(bins{j}))'];
      stairLabels{j} = sprintf('Region %d', j);
    end

    stairs([1:BIN_COUNT], stairData);
    legend(stairLabels);
    title(sprintf('Unconnected Feature %03d: %s', i, weightNames{i}));


    saveas(gcf, sprintf('feature%03d.pdf', i));
    disp(sprintf('feature%03d.pdf', i));
  end

end