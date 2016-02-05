function makeSegMovie( segmentation, raw, outputFileAdress, nrElementsInSegmentation)
% Make movie of a segmented part of the dataset

if nargin < 4
    % -1 to not count zero as an element
    nrElementsInSegmentation = length(unique(segmentation(:))) - 1;
end

% Generate distinguishable colors, avoid white and black
colors = distinguishable_colors(nrElementsInSegmentation, [0 0 0; 1 1 1]);

writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 16;
open(writerObj);
for f=1:size(raw,3)
    thisRaw = repmat(raw(:,:,f),1,1,3);
    thisSeg = label2rgb(segmentation(:,:,f), colors, [0 0 0]);
    frame = imfuse(thisRaw, thisSeg);
    writeVideo(writerObj,frame);
end
close(writerObj);

end