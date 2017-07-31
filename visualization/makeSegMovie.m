function makeSegMovie( segmentation, raw, outputFileAdress,outputFileAdress2,nrElementsInSegmentation)
% Make movie of a segmented part of the dataset

if nargin < 5
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
writerObj2 = VideoWriter(outputFileAdress2);
writerObj2.FrameRate = 10;
open(writerObj2);

for f=1:size(raw,3)
    thisRaw = repmat(raw(:,:,f),1,1,3); % single
    thisSeg = label2rgb(segmentation(:,:,f), colors, [0 0 0]); % uint8
    frame = imfuse(thisRaw, thisSeg, 'blend'); % uint8
    % rotate and flip to match wK view 
    frame2 = horzcat( flipdim(imrotate(im2uint8(thisRaw),90),1), flipdim(imrotate(frame,90),1)); % uint8
    writeVideo(writerObj,frame);
    writeVideo(writerObj2,frame2);
end
close(writerObj);
close(writerObj2);

end
