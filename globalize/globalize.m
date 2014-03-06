function [ corrs ] = globalize( cubeCoords1, cubeCoords2 )
%GLOBALIZE Summary of this function goes here
%   Detailed explanation goes here


cubeSize = 148;
overlap = 20;

seg_path = '/x/CortexConnectomics/Manuel/seg';
cube1 = readKnossosCube(seg_path,'2012-09-28_ex145_07x2_mag1',cubeCoords1 ,'uint16','', 'raw', 148);    
cube2 = readKnossosCube(seg_path,'2012-09-28_ex145_07x2_mag1',cubeCoords2 ,'uint16','', 'raw', 148);    
coordShift = cubeCoords2 - cubeCoords1;

overlap1 = cube1((1+max(0,coordShift(1))*(cubeSize-overlap)): (cubeSize+min(0,coordShift(1))*(cubeSize-overlap)),...
    (1+max(0,coordShift(2))*(cubeSize-overlap)): (cubeSize+min(0,coordShift(2))*(cubeSize-overlap)),...
    (1+max(0,coordShift(3))*(cubeSize-overlap)): (cubeSize+min(0,coordShift(3))*(cubeSize-overlap)));
overlap2 = cube2((1-min(0,coordShift(1))*(cubeSize-overlap)): (cubeSize-max(0,coordShift(1))*(cubeSize-overlap)),...
    (1-min(0,coordShift(2))*(cubeSize-overlap)): (cubeSize-max(0,coordShift(2))*(cubeSize-overlap)),...
    (1-min(0,coordShift(3))*(cubeSize-overlap)): (cubeSize-max(0,coordShift(3))*(cubeSize-overlap)));

props1 = regionprops(overlap1, 'PixelIdxList');
props2 = regionprops(overlap2, 'PixelIdxList');

objects1 = getObjects(props1,overlap1);
objects2 = getObjects(props2,overlap2);

corrs = [];
for i = 1:size(objects1,2)
    voxSize1 = size(objects1(i).PixelIdxList,1);
    for j = 1:size(objects2,2)
        voxSize2 = size(objects2(j).PixelIdxList,1);
        commonVox = size(intersect(objects1(i).PixelIdxList,objects2(j).PixelIdxList),1);
        ratio = commonVox / voxSize1;
        if(commonVox > 0 && ratio > 0.1)
            corrs = [corrs; objects1(i).col  objects2(j).col  voxSize1  voxSize2  commonVox ];
        end
    end
end

end


function objects = getObjects(props,overlap)
% discard empty objects
count = 0;
for j = 1:size(props,1)
    if(size(props(j).PixelIdxList,1) ~= 0)
        count = count+1;
        objects(count) = props(j);
    end
end

% add color values 
tmp=cell(size(objects));
[objects(:).col]=deal(tmp{:});
for j = 1:size(objects,2)
    [x, y, z] = ind2sub(size(overlap),objects(j).PixelIdxList(1,1));
    objects(j).col = overlap(x,y,z);
end
end

