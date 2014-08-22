function localVsGlobalSegmentation(p)

raw = readKnossosRoi('/zdata/manuel/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/', '2012-09-28_ex145_07x2_corrected_mag1', [641 8320; 769 5888; 1000 1000]);    
segmentation = readKnossosRoi('/zdata/manuel/results/pipeline/20140805T180102/seg/', '2012-09-28_ex145_07x2_corrected_mag1', [641 8320; 769 5888; 1000 1000], 'uint32', '', 'raw');

im = repmat(raw,[1 1 3]);
imwrite(im, '/zdata/manuel/sync/largeRaw.tif');

% Join according to correspondences
load([p.saveFolder 'seg/m.mat']);
for i=1:length(components)
    for j=2:length(components{i})
        segmentation(segmentation == components{i}(j)) = components{i}(1);
    end
end

% renumber segments
uniqueID = unique(segmentation(:));
segNew = segmentation;
for i=1:length(uniqueID)
    segNew(segmentation == uniqueID(i)) = i;
end
segNew = segNew - 1;

uniqueColors = distinguishable_colors(100, [1 1 1; 0 0 0]);
cm = [0 0 0; repmat(uniqueColors,[1000,1])];
colors = label2rgb(segNew, cm);
imwrite(colors, '/zdata/manuel/sync/largeSeg.tif');
toStore = im * 0.8 + colors * 0.2;
imwrite(toStore, '/zdata/manuel/sync/largeSegOverlay.tif');

end
