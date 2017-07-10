function skeletonFromAgglo(edges, segmentMeta, classes, name, outputFolder)

    if ~isempty(classes)
        skeletonNames = arrayfun(@(x)[name num2str(x, '%.4i')], 1:length(classes), 'uni', 0);
        connectEM.generateSkeletonFromAgglo(edges, segmentMeta.point', ...
            classes, skeletonNames, ...
            outputFolder, segmentMeta.maxSegId);
    end

end

