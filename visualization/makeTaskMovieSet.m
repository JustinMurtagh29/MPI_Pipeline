function makeTaskMovieSet(raw, seg, outputFileAdress)
% Raw with one object seg overlaid in red and saved as movie
% Without need for graphical display (e.g. cluster)

writerObj = VideoWriter(outputFileAdress);
writerObj.FrameRate = 15;
open(writerObj);
frameWithoutStart = 0;
for f=1:size(raw,3)
    thisRaw = raw(:,:,f);
    thisSeg = seg(:,:,f);
    if ~any(thisSeg(:))
        frameWithoutStart = frameWithoutStart + 1;
        if frameWithoutStart == 10
            break;
        end
    end
    thisRed = thisRaw;
    thisRed(thisSeg) = thisRed(thisSeg) + 20;
    thisRed(thisRed > 255) = 255;
    frame = zeros([size(thisRaw) 3], 'uint8');
    frame(:,:,1) = thisRed;
    frame(:,:,2) = thisRaw;
    frame(:,:,3) = thisRaw;
    writeVideo(writerObj,frame);
end
close(writerObj);

end
