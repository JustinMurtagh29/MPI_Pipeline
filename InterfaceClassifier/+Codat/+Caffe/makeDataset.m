function makeDataset( raw, target, folder, imSize, numClassInstances, augmentData, imageFormat )
%MAKEDATASET Make dataset for caffe.

fprintf('Sampling images.\n');
[ X, y, coords ] = Codat.Caffe.sampleImages( raw, target, imSize, numClassInstances, augmentData);

if ~exist(folder,'dir')
    fprintf('Creating folder %s.\n',folder);
    mkdir(folder)
end

fprintf('Saving to folder.\n');
folderContent = dir(folder);
files = {folderContent(3:end).name};
noImageFiles = sum(cellfun(@(x)~isempty(strfind(x,imageFormat)),files));
fid = fopen(sprintf('%s%simages.txt',folder,filesep),'at+');
for i = 1:length(y)
    imwrite(X(:,:,i),sprintf('%s%s%d.%s', folder, filesep, i + noImageFiles, imageFormat));
    fprintf(fid,'%d.%s %d\n', i + noImageFiles, imageFormat, y(i));
end
fclose(fid);
m = matfile([folder filesep 'imageCoords'],'Writable',true);
try %append to already existing file
    coordsOld = m.coords;
    coords = [coords,coordsOld];
catch
end
m.coords = coords;
clear
end

