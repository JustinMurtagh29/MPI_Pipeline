function localVsGlobalSegmentationvideo(p)

%raw = readKnossosRoi('/zdata/manuel/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/', '2012-09-28_ex145_07x2_corrected_mag1', [641 8320; 769 5888; 769 1024]);    
%segmentation = readKnossosRoi('/zdata/manuel/results/pipeline/20140805T180102/seg/', '2012-09-28_ex145_07x2_corrected_mag1', [641 8320; 769 5888; 769 1024], 'uint32', '', 'raw');

raw = readKnossosRoi('/zdata/manuel/data/cortex/2012-09-28_ex145_07x2_corrected/color/1/', '2012-09-28_ex145_07x2_corrected_mag1', [641 8320; 769 5888; 769 779]);    
segmentation = readKnossosRoi('/zdata/manuel/results/pipeline/20140805T180102/seg/', '2012-09-28_ex145_07x2_corrected_mag1', [641 8320; 769 5888; 769 779], 'uint32', '', 'raw');

% Join according to correspondences
load([p.saveFolder 'seg/m.mat']);
for i=1:length(components)
    tic;
    display(['Joining correspondence component ' num2str(i) ' of ' num2str(length(components))]);
    for j=2:length(components{i})
        segmentation(segmentation == components{i}(j)) = components{i}(1);
    end
    toc;
end

% renumber segments
uniqueID = unique(segmentation(:));
segNew = segmentation;
for i=1:length(uniqueID) 
    display(['Renumbering segment ' num2str(i) ' of ' num2str(length(uniqueID))]);
    segNew(segmentation == uniqueID(i)) = i;
end
segNew = segNew - 1;

% Writing to video file
uniqueColors = distinguishable_colors(100, [1 1 1; 0 0 0]);
cm = repmat(uniqueColors,[10000,1]);
writerObj = VideoWriter('/zdata/manuel/sync/largeSegMovie.avi');
writerObj.FrameRate = 10;
open(writerObj);
for i=1:size(raw,3)
    display(['Rendering frame ' num2str(i) ' of ' num2str(size(raw,3))]);    
    grey = repmat(raw(:,:,i), [1 1 3]);
    colors = label2rgb(segNew(:,:,i), cm, [0 0 0]);
    toStore = grey * 0.85 + colors * 0.15;
    writeVideo(writerObj,toStore);
end
close(writerObj);

end
