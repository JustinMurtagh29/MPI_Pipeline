function findNodes(segFile, tileBorder, lToGFile, segmentFile)

    % Load local segmentation
    [~,segSmall] = loadSegData(segFile, tileBorder);
    % Load local to global mapping
    %load(lToGFile);
    
    % Extract linear idx for each segment and mean intensity
    segments = regionprops(segSmall, 'PixelIdxList');
    for i=1:length(segments)
        segments(i).Id = i;
        %    segments(i).Id = globalIds(i == localIds);
    end
    segments(arrayfun(@(x)isempty(x.PixelIdxList), segments)) = [];
    save(segmentFile, 'segments');

end
