function makeSegMovie( segmentation, raw, outputFileAdress, nrElementsInSegmentation)
% Make movie of a segmented part of the dataset

if nargin < 4
    % -1 to not count zero as an element
    nrElementsInSegmentation = length(unique(segmentation(:))) - 1;
end

% Restrict to 500 colors if more segments present (hard to distinguish more)
if nrElementsInSegmentation > 500
    nrColors = 500;
else
    nrColors = nrElementsInSegmentation;
end

% Generate distinguishable colors, avoid white and black
colors = distinguishable_colors(nrColors, [0 0 0; 1 1 1]);
colors = repmat(colors, ceil(nrElementsInSegmentation/nrColors), 1);

writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 10;
open(writerObj);
for f=1:size(raw,3)
    thisRaw = repmat(raw(:,:,f),1,1,3);
    thisSeg = label2rgb(segmentation(:,:,f), colors, [0 0 0]);
    frame = imfuse(thisRaw, thisSeg, 'blend');
    writeVideo(writerObj,frame);
end
close(writerObj);

end
