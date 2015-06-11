FileTif='/home/mberning/fsHest/Data/huay/for MB/m150416_013.tif';
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
 
FinalImage=zeros(nImage,mImage,NumberImages,'uint8');
for i=1:NumberImages
   FinalImage(:,:,i)=imread(FileTif,'Index',i);
end

%% Seperate channels
greenChannel = FinalImage(:,:,1:2:end-1);
redChannel = FinalImage(:,:,2:2:end);
for i=1:size(greenChannel,3)
    combinedChannels(:,:,1) = redChannel(:,:,i) ./ 255;
    combinedChannels(:,:,2) = greenChannel(:,:,i) ./ 255;
    combinedChannels(:,:,3) = 0;
    imwrite(combinedChannels, ['/home/mberning/Desktop/YHtest/' num2str(i) '.tif']);
end


