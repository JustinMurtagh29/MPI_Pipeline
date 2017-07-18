function skeletonFromAgglo(edges, segmentMeta, classes, name, outputFolder)
    if ~exist(outputFolder,'dir')
	mkdir(outputFolder);
    end
    % What did you do Alessandro :D
    if size(segmentMeta.point,2) ~= 3
        point = segmentMeta.point';
    else
        point = segmentMeta.point;
    end

    if ~isempty(classes)
        skeletonNames = arrayfun(@(x)[name num2str(x, '%.4i')], 1:length(classes), 'uni', 0);
        connectEM.generateSkeletonFromAgglo(edges, point, ...
            classes, skeletonNames, ...
            outputFolder, segmentMeta.maxSegId);
    end

end

