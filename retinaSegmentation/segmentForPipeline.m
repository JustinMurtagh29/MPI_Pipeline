function segmentForPipelineRetina( root,folders, prefix, bbox, segFunction, saveFile )



% Load classification from three affines fodlers
affx = loadClassData([root folders{1}],prefix,bbox);
affy = loadClassData([root folders{2}],prefix,bbox);
affz = loadClassData([root folders{3}],prefix,bbox);

%morphR done here
aff = {affx affy affz};
clear affx affy affz

rT=1;
[x,y,z] = meshgrid(-rT:rT,-rT:rT,-rT:rT);
se = (x/rT).^2 + (y/rT).^2 + (z/rT).^2 <= 1;

for dir=1:length(aff)
     
    % Opening by reconstruction
    affRecon = imerode(aff{dir}, se);
    affRecon = imreconstruct(affRecon, aff{dir});
    affRecon = imcomplement(affRecon);

    % Closing by reconstruction
    affReconDilated = imdilate(affRecon, se);
    affReconDilated = imcomplement(affReconDilated);

    aff{dir} = imreconstruct(affReconDilated, affRecon);
end

% Perform segmentation
seg = segFunction(aff);
seg = uint16(seg{1,1});

clear affReconRecon
% If folder does not exist, create it
saveFolder = fileparts(saveFile);
if ~exist(saveFolder, 'dir')
    mkdir(saveFolder);
end

% Save segmentation to MATLAB file in 'saveFolder'
save(saveFile, 'seg');

end
