function makeLCMovie( path )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

files = dir([path '*.tif']);
writerObj = VideoWriter('C:\Users\mberning\Desktop\taskMovie.avi');
writerObj.FrameRate = 16;
open(writerObj);
for i=1:length(files)
        writeVideo(writerObj,imread([path files(i).name]));
end
close(writerObj);
close all;

end