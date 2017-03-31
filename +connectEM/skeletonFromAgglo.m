function skeletonFromAgglo(edges, segmentMeta, classes, name, outputFolder)

    skeletonNames = cellfun(@(x)[name num2str(x, '%5i')], 1:length(classes), 'uni', 0);
    connectEM.generateSkeletonFromAgglo(edges, segmentMeta.point, ...
        classes, skeletonNames, ...
        outputFolder, segmentMeta.maxSegId);

end

